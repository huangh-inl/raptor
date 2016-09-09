#include "ThermiteMaterial.h"

template<>
InputParameters validParams<ThermiteMaterial>()
{
  InputParameters params = validParams<Material>();

  params.addCoupledVar("temp_var", "Temperature at current point");

  params.addParam<std::vector<std::string> >("reaction_system", "names of reactants and products");
  params.addParam<std::vector<Real> >("sto_coeff", "stoichiometry coefficients for reactans and products");
  params.addParam<std::vector<Real> >("molecular_weights", "molecule weights for reactants and products kg/mol");

  params.addParam<std::vector<Real> >("melting_temp", "metlting temperature for reactants and products in Kelvin");
  params.addParam<std::vector<Real> >("vaporization_temp", "boiling temperature for reactants and products in Kelvin");
  
  params.addParam<std::vector<Real> >("latent_s_to_l", "latent heat from soild to liquid in KJ/Kg");
  params.addParam<std::vector<Real> >("latent_l_to_g", "latent heat from liquid to gas in KJ/Kg");

  params.addParam<Real>("ignition_temp", 1100, "ignition temperature in Kelvin");
  params.addParam<Real>("rate_constant", 1.0,  "reaction rate");

  
  params.addParam<Real>("initial_thermal_conductivity", 1.0, "The diffusion coefficient of the Reactants kW/(m.K) ");
  params.addParam<Real>("initial_specific_heat", 1.0, "The specific heat of the Reactants KJ/(kg.K)");
  params.addParam<Real>("initial_density", 1.0, "The density the Reactants kg/m^3");

  return params;
}

ThermiteMaterial::ThermiteMaterial(const InputParameters & parameters) :
  Material(parameters),

  // Get a parameter value for the specific heat
  _initial_thermal_conductivity(getParam<Real>("initial_thermal_conductivity")),
  _initial_specific_heat(getParam<Real>("initial_specific_heat")),
  _initial_density(getParam<Real>("initial_density")),

  _ignition_temp(getParam<Real>("ignition_temp")),
  _kr(getParam<Real>("rate_constant")),

  // Declare that this material is going to have a real
  // value property named "specific_heat" that Kernels can use.
  _thermal_conductivity(declareProperty<Real>("thermal_conductivity")),  
  _specific_heat(declareProperty<Real>("specific_heat")),  
  _density(declareProperty<Real>("density")),  

  _r_ex(declareProperty<Real>("reaction_extent")),
  _r_ex_old(declarePropertyOld<Real>("reaction_extent")),
  _temp(coupledValue("temp_var")) 
{
  _species = getParam<std::vector<std::string>>("reaction_system");
  _sto_coeff = getParam<std::vector<Real>>("sto_coeff");
  _molecular_weights = getParam<std::vector<Real>>("molecular_weights");
  _melting_temp = getParam<std::vector<Real>>("melting_temp");
  _vaporization_temp = getParam<std::vector<Real>>("vaporization_temp");
  _latent_s_to_l = getParam<std::vector<Real>>("latent_s_to_l");
  _latent_l_to_g = getParam<std::vector<Real>>("latent_l_to_g");

  unsigned int n(_species.size());
  _initial_mass_fraction.resize(n);
  _current_mass_fraction.resize(n);
  _initial_mass_conc.resize(n);
  _current_mass_conc.resize(n);
  _initial_molar_conc.resize(n);
  _current_molar_conc.resize(n);
  
  //right now we assume reactants are stoichiometrically matched for complete reaction
  for (unsigned int i = 0; i < _species.size(); ++i)
  {
    if ( i <= 1)
    {
      _initial_mass_fraction[i] = _sto_coeff[i] * _molecular_weights[i] /
                                ( _sto_coeff[0] *  _molecular_weights[0] +  _sto_coeff[1] *  _molecular_weights[1]);
      _initial_mass_conc[i] = _initial_density * _initial_mass_fraction[i];
    }
    else
    {
      _initial_mass_fraction[i] = 0.0; //no reaction products at the begining of simulation
      _initial_mass_conc[i] = 0.0;
    }
    
    _current_mass_fraction[i] = _initial_mass_fraction[i];
    _current_mass_conc[i] = _initial_mass_conc[i];

    _initial_molar_conc[i] = _initial_mass_conc[i] / _molecular_weights[i];
    _current_molar_conc[i] = _initial_molar_conc[i];
  }
}

void
ThermiteMaterial::initQpStatefulProperties()
{
  _r_ex[_qp] = 0.0;
  _r_ex_old[_qp] = 0.0;
}


void
ThermiteMaterial::computeQpProperties()
{

  Real T = _temp[_qp];
  Real x_old = _r_ex_old[_qp];
  Real dx(0.0);
  
  if ( T < _ignition_temp)
    dx = 0.0;
  else if (T >= _ignition_temp && x_old < 1.0)
    dx = _kr * _dt / (_kr * _dt + 1.0)  * (1.0 - x_old);
  else if (x_old >= 1.0)
    dx = 0.0;
  else
    dx = 0.0;

//  if (dx > 0) std::cout<<dx<<std::endl;

  Real x = x_old + dx;
  if (x < 0.0)
    x = 0.0;
  if (x > 1.0)
    x = 1.0;

  
  _r_ex[_qp] = x;

//  if (_r_ex[_qp]  > 0) std::cout<< _r_ex[_qp] <<std::endl;

// current mass fraction of species
  //reacted oxide in mole/m^3
  Real reacted_amount =  _initial_molar_conc[0] * x;
  Real incre_reacted_amount = _initial_molar_conc[0] * dx;
  
  Real total_mass_conc(0.0);
  
  for (unsigned int i = 0; i < _species.size(); ++i)
  {
    if (i <= 1) //reactants
    {
      _current_molar_conc[i] = _initial_molar_conc[i] - _sto_coeff[i] / _sto_coeff[0] * reacted_amount;
    }
    
    else //products
    {
      _current_molar_conc[i] = _sto_coeff[i] / _sto_coeff[0] * reacted_amount;
    }
    
    _current_mass_conc[i] =  _current_molar_conc[i] * _molecular_weights[i];
    total_mass_conc += _current_mass_conc[i];
  }
  
  // for now the thermal conductivity and density are treated as constants
  _thermal_conductivity[_qp] =  _initial_thermal_conductivity; 
  _density[_qp] = _initial_density;

  // loop over all species, check phases
  std::vector<Real> Cp;
  Cp.resize(_species.size());
    
  for (unsigned int i = 0; i < _species.size(); ++i)
  {
    Cp[i] = _initial_specific_heat;
    if ( T >= _melting_temp[i] && T < _vaporization_temp[i])
      Cp[i] += _latent_s_to_l[i];
    if ( T >=  _vaporization_temp[i])
      Cp[i] += _latent_l_to_g[i];
  }
  
  Real effective_Cp(0.0);
    
  for (unsigned int i = 0; i < _species.size(); ++i)
    effective_Cp += Cp[i] * _current_mass_conc[i] / total_mass_conc;
    
  _specific_heat[_qp] = effective_Cp;
    
  //  std::cout<<_thermal_conductivity[_qp]<<std::endl;
}
