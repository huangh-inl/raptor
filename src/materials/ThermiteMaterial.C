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
  params.addParam<std::vector<Real> >("specific_heat_solid", "specific heat of all component in solid KJ/(Kg.K)");
  params.addParam<std::vector<Real> >("specific_heat_liquid", "specific heat of all component in liquid KJ/(Kg.K)");
  params.addParam<std::vector<Real> >("component_density","density of all component Kg/m^3");
  params.addParam<std::vector<Real> >("component_thermal_cond","thermal conductivity of all component kW/(m.K)");

  params.addParam<Real>("ignition_temp", 1100, "ignition temperature in Kelvin");
  params.addParam<Real>("rate_constant", 1.0,  "reaction rate");

  
  params.addParam<Real>("initial_thermal_conductivity", 1.0, "The diffusion coefficient of the Reactants kW/(m.K) ");
  params.addParam<Real>("initial_specific_heat", 1.0, "The specific heat of the Reactants KJ/(kg.K)");
  params.addParam<Real>("initial_density", 1.0, "The density the Reactants kg/m^3");
  params.addParam<Real>("heat_fraction", 1.0, "The heat fraction used in melting and heat diffusion, 0-1");

  return params;
}

ThermiteMaterial::ThermiteMaterial(const InputParameters & parameters) :
  Material(parameters),

  // Get a parameter value for the specific heat
  _initial_thermal_conductivity(getParam<Real>("initial_thermal_conductivity")),
  _initial_specific_heat(getParam<Real>("initial_specific_heat")),
  _initial_density(getParam<Real>("initial_density")),
  _heat_fraction(getParam<Real>("heat_fraction")),

  _ignition_temp(getParam<Real>("ignition_temp")),
  _kr(getParam<Real>("rate_constant")),

  // Declare that this material is going to have a real
  // value property named "specific_heat" that Kernels can use.
  _thermal_conductivity(declareProperty<Real>("thermal_conductivity")),  
  _specific_heat(declareProperty<Real>("specific_heat")),  
  _density(declareProperty<Real>("density")),  

  _r_ex(declareProperty<Real>("reaction_extent")),
  _r_ex_old(declarePropertyOld<Real>("reaction_extent")),

  _heat_sink_melting(declareProperty<Real>("melting_heatsink")),
  _heat_sink_boiling(declareProperty<Real>("boiling_heatsink")),

  _react1_status(declareProperty<Real>("Al_phase")),
  _react2_status(declareProperty<Real>("Fe2O3_phase")),
  _prod1_status(declareProperty<Real>("Al2O3_phase")),
  _prod2_status(declareProperty<Real>("Fe_phase")),

  _react1_phase(declareProperty<Real>("Al_state")),
  _react2_phase(declareProperty<Real>("Fe2O3_state")),
  _react3_phase(declareProperty<Real>("Al2O3_state")),
  _react4_phase(declareProperty<Real>("Fe_state")),

  _react1_solid(declareProperty<Real>("Al_solid")),
  _react1_solid_old(declarePropertyOld<Real>("Al_solid")),
  _react1_liquid(declareProperty<Real>("Al_liquid")),
  _react1_liquid_old(declarePropertyOld<Real>("Al_liquid")),
  _react1_gas(declareProperty<Real>("Al_gas")),
  _react1_gas_old(declarePropertyOld<Real>("Al_gas")),

  _react2_solid(declareProperty<Real>("Fe2O3_solid")),
  _react2_solid_old(declarePropertyOld<Real>("Fe2O3_solid")),
  _react2_liquid(declareProperty<Real>("Fe2O3_liquid")),
  _react2_liquid_old(declarePropertyOld<Real>("Fe2O3_liquid")),
  _react2_gas(declareProperty<Real>("Fe2O3_gas")),
  _react2_gas_old(declarePropertyOld<Real>("Fe2O3_gas")),

  _react3_solid(declareProperty<Real>("Al2O3_solid")),
  _react3_solid_old(declarePropertyOld<Real>("Al2O3_solid")),
  _react3_liquid(declareProperty<Real>("Al2O3_liquid")),
  _react3_liquid_old(declarePropertyOld<Real>("Al2O3_liquid")),
  _react3_gas(declareProperty<Real>("Al2O3_gas")),
  _react3_gas_old(declarePropertyOld<Real>("Al2O3_gas")),

  _react4_solid(declareProperty<Real>("Fe_solid")),
  _react4_solid_old(declarePropertyOld<Real>("Fe_solid")),
  _react4_liquid(declareProperty<Real>("Fe_liquid")),
  _react4_liquid_old(declarePropertyOld<Real>("Fe_liquid")),
  _react4_gas(declareProperty<Real>("Fe_gas")),
  _react4_gas_old(declarePropertyOld<Real>("Fe_gas")),
  
  _temp(coupledValue("temp_var")) 
{
  _species = getParam<std::vector<std::string>>("reaction_system");
  _sto_coeff = getParam<std::vector<Real>>("sto_coeff");
  _molecular_weights = getParam<std::vector<Real>>("molecular_weights");
  _melting_temp = getParam<std::vector<Real>>("melting_temp");
  _vaporization_temp = getParam<std::vector<Real>>("vaporization_temp");
  _latent_s_to_l = getParam<std::vector<Real>>("latent_s_to_l");
  _latent_l_to_g = getParam<std::vector<Real>>("latent_l_to_g");
  _specific_heat_solid = getParam<std::vector<Real>>("specific_heat_solid");
  _specific_heat_liquid = getParam<std::vector<Real>>("specific_heat_liquid");
  _component_density = getParam<std::vector<Real>>("component_density");
  _component_thermal_cond = getParam<std::vector<Real>>("component_thermal_cond");
  

  unsigned int n(_species.size());
  _initial_mass_fraction.resize(n);
  _current_mass_fraction.resize(n);
  _current_vol_fraction.resize(n);
  _initial_mass_conc.resize(n);
  _current_mass_conc.resize(n);
  _initial_molar_conc.resize(n);
  _current_molar_conc.resize(n);
  
  // Those parameters are calculated from the comparison of released heat from reaction and latent heat
  _reacted_mass_conc.resize(n);

  Real tot_mass(0.0);
  
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

    tot_mass += _current_mass_conc[i];  // Total mass concentration

    _initial_molar_conc[i] = _initial_mass_conc[i] / _molecular_weights[i];
    _current_molar_conc[i] = _initial_molar_conc[i];
  }  

  for (unsigned int i = 0; i < _species.size(); ++i)
  {
    _current_vol_fraction[i] = (_current_mass_conc[i]/_component_density[i])/(tot_mass/_initial_density);
  }
  
}

void
ThermiteMaterial::initQpStatefulProperties()
{
  _r_ex[_qp] = 0.0;
  _r_ex_old[_qp] = 0.0;
  
  _heat_sink_melting[_qp] = 0.0;
  _heat_sink_boiling[_qp] = 0.0;

  _react1_solid[_qp] = _initial_mass_conc[0];
  _react1_solid_old[_qp] = _initial_mass_conc[0];
  _react1_liquid[_qp] = 0.0;
  _react1_liquid_old[_qp] = 0.0;
  _react1_gas[_qp] = 0.0;
  _react1_gas_old[_qp] = 0.0;

  _react2_solid[_qp] = _initial_mass_conc[1];
  _react2_solid_old[_qp] = _initial_mass_conc[1];
  _react2_liquid[_qp] = 0.0;
  _react2_liquid_old[_qp] = 0.0;
  _react2_gas[_qp] = 0.0;
  _react2_gas_old[_qp] = 0.0;

  _react3_solid[_qp] = _initial_mass_conc[2];
  _react3_solid_old[_qp] = _initial_mass_conc[2];
  _react3_liquid[_qp] = 0.0;
  _react3_liquid_old[_qp] = 0.0;
  _react3_gas[_qp] = 0.0;
  _react3_gas_old[_qp] = 0.0;
  
  _react4_solid[_qp] = _initial_mass_conc[3];
  _react4_solid_old[_qp] = _initial_mass_conc[3];
  _react4_liquid[_qp] = 0.0;
  _react4_liquid_old[_qp] = 0.0;
  _react4_gas[_qp] = 0.0;
  _react4_gas_old[_qp] = 0.0;
  
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
  {
    x = 0.0;
  }
  if (x > 1.0)
  {  
    x = 1.0;
  }
  
  _r_ex[_qp] = x;
  dx = _r_ex[_qp] - _r_ex_old[_qp];
  
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

  // Calculate the reacted mass amount and produced mass amount at each time
  for (unsigned int i = 0; i < _species.size(); ++i)
  {
    _reacted_mass_conc[i] = (_sto_coeff[i] / _sto_coeff[0]) * incre_reacted_amount * _molecular_weights[i];
  }
  

  // Update the mixture density with the component mass change
  Real rho(0.0);
    
  for (unsigned int i = 0; i < _species.size(); ++i)
  {
    rho += (_current_mass_conc[i] / total_mass_conc)/_component_density[i];
  }

  
 //********************************************************************
 // Update all the solid,liquid,gas phase conditions due to reaction
 //********************************************************************

  // Reactant 1
  Real r1_solid_old(0.0);
  Real r1_solid(0.0);
  Real r1_liquid_old(0.0);
  Real r1_liquid(0.0);
  Real r1_gas_old(0.0);
  Real r1_gas(0.0);
  Real r1_mass(0.0);
  
  r1_solid_old = _react1_solid_old[_qp];
  r1_liquid_old = _react1_liquid_old[_qp];
  r1_gas_old = _react1_gas_old[_qp];
  r1_mass = _reacted_mass_conc[0];

  if (r1_solid_old > 0)
  {
    if (r1_solid_old >= r1_mass)
    {  
      r1_solid = r1_solid_old - r1_mass;
      r1_liquid = r1_liquid_old;
      r1_gas = r1_gas_old;
    }
    
    else
    {
      r1_liquid = r1_liquid_old - (r1_mass - r1_solid_old);
      r1_solid = 0.0;
      r1_gas = r1_gas_old;
    }   
  }
  
  else if (r1_solid_old == 0 && r1_liquid_old > 0)
  {
    if (r1_liquid_old >= r1_mass)
    {
      r1_liquid = r1_liquid_old - r1_mass;
      r1_solid = 0.0;
      r1_gas = r1_gas_old;
    }
    else
    {
      r1_gas = r1_gas_old - (r1_mass - r1_liquid_old);
      r1_liquid = 0.0;
      r1_solid = 0.0;
    }
    
  }

  else if (r1_solid_old == 0 && r1_liquid_old == 0)
  {
    r1_gas = r1_gas_old - r1_mass;
    r1_solid = 0.0;
    r1_liquid = 0.0;
  }  

  // Reactant 2
  Real r2_solid_old(0.0);
  Real r2_solid(0.0);
  Real r2_liquid_old(0.0);
  Real r2_liquid(0.0);
  Real r2_gas_old(0.0);
  Real r2_gas(0.0);
  Real r2_mass(0.0);
  
  r2_solid_old = _react2_solid_old[_qp];
  r2_liquid_old = _react2_liquid_old[_qp];
  r2_gas_old = _react2_gas_old[_qp];
  r2_mass = _reacted_mass_conc[1];

  if (r2_solid_old > 0)
  {
    if (r2_solid_old >= r2_mass)
    {
      r2_solid = r2_solid_old - r2_mass;
      r2_liquid = r2_liquid_old;
      r2_gas = r2_gas_old;
    }
    else
    {
      r2_liquid = r2_liquid_old - (r2_mass - r2_solid_old);
      r2_solid = 0.0;
      r2_gas = r2_gas_old;
    }   
  }
  
  else if (r2_solid_old == 0 && r2_liquid_old > 0)
  {
    if (r2_liquid_old >= r2_mass)
    {
      r2_liquid = r2_liquid_old - r2_mass;
      r2_solid = 0.0;
      r2_gas = r2_gas_old;
    }
    else
    {
      r2_gas = r2_gas_old - (r2_mass - r2_liquid_old);
      r2_liquid = 0.0;
      r2_solid = 0.0;
    }
    
  }

  else if (r2_solid_old == 0 && r2_liquid_old == 0)
  {  
    r2_gas = r2_gas_old - r2_mass;
    r2_solid = 0.0;
    r2_liquid = 0.0;
  }
  

  //Reactant 3 - Products
  Real r3_solid_old(0.0);
  Real r3_solid(0.0);
  Real r3_liquid_old(0.0);
  Real r3_liquid(0.0);
  Real r3_gas_old(0.0);
  Real r3_gas(0.0);
  Real r3_mass(0.0);  
  
  r3_solid_old = _react3_solid_old[_qp];
  r3_liquid_old = _react3_liquid_old[_qp];
  r3_gas_old = _react3_gas_old[_qp];
  r3_mass = _reacted_mass_conc[2];

  r3_solid = r3_solid_old + r3_mass;
  r3_liquid = r3_liquid_old;
  r3_gas = r3_gas_old;

  // Reactant 4 - Products
  Real r4_solid_old(0.0);
  Real r4_solid(0.0);
  Real r4_liquid_old(0.0);
  Real r4_liquid(0.0);
  Real r4_gas_old(0.0);
  Real r4_gas(0.0);
  Real r4_mass(0.0);  
  
  r4_solid_old = _react4_solid_old[_qp];
  r4_liquid_old = _react4_liquid_old[_qp];
  r4_gas_old = _react4_gas_old[_qp];
  r4_mass = _reacted_mass_conc[3];

  r4_solid = r4_solid_old + r4_mass;
  r4_liquid = r4_liquid_old;
  r4_gas = r4_gas_old;

  // Calculate the volume fraction of each component
  for (unsigned int i = 0; i < _species.size(); ++i)
  {
    _current_vol_fraction[i] = (_current_mass_conc[i]/_component_density[i])/(total_mass_conc/(1/rho));
  }
  
// // // ***********************************************************
// // //   Update the heat sink term according to the reaction
// // // ***********************************************************

  Real current_r_heat = _heat_fraction * 5322.746 * 4360 * (_r_ex[_qp]-_r_ex_old[_qp])/_dt;  // current released reaction heat
  Real ini_r_heat = current_r_heat;

  Real r1_melting_req(0.0);  // melting heat - required heat for totally melting
  Real r1_melting_eff(0.0);  // melting heat - effective heat
  Real r1_boiling_req(0.0);  // boiling heat - required
  Real r1_boiling_eff(0.0);  // boiling heat
  Real r1_solid_new(0.0);
  Real r1_liquid_new(0.0);
  Real r1_gas_new(0.0);
  Real r1_melting_conc(0.0);
  
  Real r2_melting_req(0.0);
  Real r2_melting_eff(0.0);
  Real r2_boiling_req(0.0);
  Real r2_boiling_eff(0.0);
  Real r2_solid_new(0.0);
  Real r2_liquid_new(0.0);
  Real r2_gas_new(0.0);
  Real r2_melting_conc(0.0);

  Real r3_melting_req(0.0);
  Real r3_melting_eff(0.0);
  Real r3_boiling_req(0.0);
  Real r3_boiling_eff(0.0);
  Real r3_solid_new(0.0);
  Real r3_liquid_new(0.0);
  Real r3_gas_new(0.0);
  Real r3_melting_conc(0.0);

  Real r4_melting_req(0.0);
  Real r4_melting_eff(0.0);
  Real r4_boiling_req(0.0);
  Real r4_boiling_eff(0.0);
  Real r4_solid_new(0.0);
  Real r4_liquid_new(0.0);
  Real r4_gas_new(0.0);
  Real r4_melting_conc(0.0);
  Real left_r_heat(0.0);
  Real r4_boiling_req2(0.0);
  Real r4_boiling_eff2(0.0);
  Real r4_boiling_conc(0.0);
  

  // For the first reactant
  
  if (T >= _melting_temp[0] && T < _vaporization_temp[0])
  {
    r1_melting_req = _latent_s_to_l[0] * r1_solid;
    if (current_r_heat >= r1_melting_req)
    {
      r1_melting_eff = r1_melting_req;
      r1_liquid_new = r1_solid + r1_liquid;
      r1_solid_new = 0.0;
      r1_gas_new = r1_gas;
    }
    else
    {
      r1_melting_eff = current_r_heat;
      r1_melting_conc = r1_melting_eff/_latent_s_to_l[0];
      r1_solid_new = r1_solid - r1_melting_conc;
      r1_liquid_new = r1_liquid + r1_melting_conc;
      r1_gas_new = r1_gas;
    }    
  }
  else if (T >= _vaporization_temp[0])
  {
    if (r1_solid > 0)
    {
      r1_boiling_req = _latent_s_to_l[0] * r1_solid;
      if (current_r_heat >= r1_boiling_req)
      {
        r1_boiling_eff = r1_boiling_req;
        r1_liquid_new = r1_solid + r1_liquid;
        r1_solid_new = 0.0;
        r1_gas_new = r1_gas;
      }
      else
      {
        r1_boiling_eff = current_r_heat;
        r1_melting_conc = r1_boiling_eff/_latent_s_to_l[0];
        r1_solid_new = r1_solid - r1_melting_conc;
        r1_liquid_new = r1_liquid + r1_melting_conc;
        r1_gas_new = r1_gas;
      }      
    }
    else if (r1_solid <= 0)
    {
      r1_boiling_req = _latent_l_to_g[0] * r1_liquid;
      if (current_r_heat >= r1_boiling_req)
      {
        r1_boiling_eff = r1_boiling_req;
        r1_gas_new = r1_gas + r1_liquid;
        r1_liquid_new = 0.0;
        r1_solid_new = 0.0;
      }
      else
      {
        r1_boiling_eff = current_r_heat;
        r1_melting_conc = r1_boiling_eff/_latent_l_to_g[0];
        r1_gas_new = r1_gas + r1_melting_conc;
        r1_liquid_new = r1_liquid - r1_melting_conc;
        r1_solid_new = 0.0;
      }
      
    }    
  }
  else
  {
    r1_solid_new = r1_solid;
    r1_liquid_new = r1_liquid;
    r1_gas_new = r1_gas;
  }

  // Update the reaction heat
  current_r_heat = current_r_heat - r1_melting_eff - r1_boiling_eff;

  //********************************************************
  // For the second reactant
  //********************************************************
  
  if (T >= _melting_temp[1] && T < _vaporization_temp[1])
  {
    r2_melting_req = _latent_s_to_l[1] * r2_solid;
    if (current_r_heat >= r2_melting_req)
    {
      r2_melting_eff = r2_melting_req;
      r2_liquid_new = r2_liquid + r2_solid;
      r2_solid_new = 0.0;
      r2_gas_new = r2_gas;
    }
    else
    {
      r2_melting_eff = current_r_heat;
      r2_melting_conc = r2_melting_eff/_latent_s_to_l[1];
      r2_solid_new = r2_solid - r2_melting_conc;
      r2_liquid_new = r2_liquid + r2_melting_conc;
      r2_gas_new = r2_gas;
    }    
  }
  else if (T >= _vaporization_temp[1])
  {
    if (r2_solid > 0)
    {
      r2_boiling_req = _latent_s_to_l[1] * r2_solid;
      if (current_r_heat >= r2_boiling_req)
      {
        r2_boiling_eff = r2_boiling_req;
        r2_liquid_new = r2_solid + r2_liquid;
        r2_solid_new = 0.0;
        r2_gas_new = r2_gas;
      }
      else
      {
        r2_boiling_eff = current_r_heat;
        r2_melting_conc = r2_boiling_eff/_latent_s_to_l[1];
        r2_solid_new = r2_solid - r2_melting_conc;
        r2_liquid_new = r2_liquid + r2_melting_conc;
        r2_gas_new = r2_gas;
      }      
    }
    else if (r2_solid <= 0)
    {
      r2_boiling_req = _latent_l_to_g[1] * r2_liquid;
      if (current_r_heat >= r2_boiling_req)
      {
        r2_boiling_eff = r2_boiling_req;
        r2_gas_new = r2_gas + r2_liquid;
        r2_liquid_new = 0.0;
        r2_solid_new = 0.0;
      }
      else
      {
        r2_boiling_eff = current_r_heat;
        r2_melting_conc = r2_boiling_eff/_latent_l_to_g[1];
        r2_gas_new = r2_gas + r2_melting_conc;
        r2_liquid_new = r2_liquid - r2_melting_conc;
        r2_solid_new = 0.0;
      }
      
    }    
  }
  else
  {
    r2_solid_new = r2_solid;
    r2_liquid_new = r2_liquid;
    r2_gas_new = r2_gas;
  }

  // Update the reaction heat
  current_r_heat = current_r_heat - r2_melting_eff - r2_boiling_eff;

  // Reactant4 - Product 2
  if (T >= _melting_temp[3] && T < _vaporization_temp[3])
  {
    r4_melting_req = _latent_s_to_l[3] * r4_solid;
    if (current_r_heat >= r4_melting_req)
    {
      r4_melting_eff = r4_melting_req;
      r4_liquid_new = r4_liquid + r4_solid;
      r4_solid_new = 0.0;
      r4_gas_new = r4_gas;
    }
    else
    {
      r4_melting_eff = current_r_heat;
      r4_melting_conc = r4_melting_eff/_latent_s_to_l[3];
      r4_solid_new = r4_solid - r4_melting_conc;
      r4_liquid_new = r4_liquid + r4_melting_conc;
      r4_gas_new = r4_gas;
    }    
  }
  else if (T >= _vaporization_temp[3])
  {
    if (r4_solid > 0)
    { 
      r4_boiling_req = _latent_s_to_l[3] * r4_solid;
      if (current_r_heat >= r4_boiling_req)
      {
        r4_boiling_eff = r4_boiling_req;
        r4_liquid_new = r4_solid + r4_liquid;
        r4_solid_new = 0.0;
        r4_gas_new = r4_gas;
        left_r_heat = current_r_heat - r4_boiling_eff;
        
        if (left_r_heat > 0) // keep melting the Fe from liquid to vapor
        {
          
          r4_boiling_req2 = _latent_l_to_g[3] * r4_liquid_new;
          if (left_r_heat >= r4_boiling_req2)
          {
//            std::cout << "case 1 " << left_r_heat << "  " << r4_boiling_req2 << "  " << r4_liquid_new << std::endl;            
            r4_boiling_eff2 = r4_boiling_req2;
            r4_gas_new += r4_liquid_new;
            r4_liquid_new = 0.0;
            r4_solid_new = 0.0;
          }
          else
          {
//            std::cout << "case 2 " << left_r_heat << "  " << r4_boiling_req2 << "  " << r4_liquid_new << std::endl;
            r4_boiling_eff2 = left_r_heat;
            r4_boiling_conc = r4_boiling_eff2/_latent_l_to_g[3];
            r4_gas_new += r4_boiling_conc;
            r4_liquid_new -= r4_boiling_conc;
            r4_solid_new = 0.0;
            
          }
          
        }
        
      }
      else
      {
        r4_boiling_eff = current_r_heat;
        r4_melting_conc = r4_boiling_eff/_latent_s_to_l[3];
        r4_solid_new = r4_solid - r4_melting_conc;
        r4_liquid_new = r4_liquid + r4_melting_conc;
        r4_gas_new = r4_gas;
      }      
    }
    else if (r4_solid == 0)
    {
      
      r4_boiling_req = _latent_l_to_g[3] * r4_liquid;
      if (current_r_heat >= r4_boiling_req)
      {
        r4_boiling_eff = r4_boiling_req;
        r4_gas_new = r4_gas + r4_liquid;
        r4_liquid_new = 0.0;
        r4_solid_new = 0.0;
      }
      else
      {
        r4_boiling_eff = current_r_heat;
        r4_melting_conc = r4_boiling_eff/_latent_l_to_g[3];
        r4_gas_new = r4_gas + r4_melting_conc;
        r4_liquid_new = r4_liquid - r4_melting_conc;
        r4_solid_new = 0.0;
      }     
    }    
  }  
  else
  {
    r4_solid_new = r4_solid;
    r4_liquid_new = r4_liquid;
    r4_gas_new = r4_gas;
  }

  // Update the reaction heat
  current_r_heat = current_r_heat - r4_melting_eff - r4_boiling_eff - r4_boiling_eff2; 
   
  // Reactant 3 - Product 1

  if (T >= _melting_temp[2] && T < _vaporization_temp[2])
  {
    r3_melting_req = _latent_s_to_l[2] * r3_solid;
    if (current_r_heat >= r3_melting_req)
    {
      r3_melting_eff = r3_melting_req;
      r3_liquid_new = r3_liquid + r3_solid;
      r3_solid_new = 0.0;
      r3_gas_new = r3_gas;
    }
    else
    {
      r3_melting_eff = current_r_heat;
      r3_melting_conc = r3_melting_eff/_latent_s_to_l[2];
      r3_solid_new = r3_solid - r3_melting_conc;
      r3_liquid_new = r3_liquid + r3_melting_conc;
      r3_gas_new = r3_gas;
    }    
  }
  else if (T >= _vaporization_temp[2])
  {
    if (r3_solid > 0)
    {
      r3_boiling_req = _latent_s_to_l[2] * r3_solid;
      if (current_r_heat >= r3_boiling_req)
      {
        r3_boiling_eff = r3_boiling_req;
        r3_liquid_new = r3_solid + r3_liquid;
        r3_solid_new = 0.0;
        r3_gas_new = r3_gas;
      }
      else
      {
        r3_boiling_eff = current_r_heat;
        r3_melting_conc = r3_boiling_eff/_latent_s_to_l[2];
        r3_solid_new = r3_solid - r3_melting_conc;
        r3_liquid_new = r3_liquid + r3_melting_conc;
        r3_gas_new = r3_gas;
      }      
    }
    else if (r3_solid <= 0)
    {
      r3_boiling_req = _latent_l_to_g[2] * r3_liquid;
      if (current_r_heat >= r3_boiling_req)
      {
        r3_boiling_eff = r3_boiling_req;
        r3_gas_new = r3_gas + r3_liquid;
        r3_liquid_new = 0.0;
        r3_solid_new = 0.0;
      }
      else
      {
        r3_boiling_eff = current_r_heat;
        r3_melting_conc = r3_boiling_eff/_latent_l_to_g[2];
        r3_gas_new = r3_gas + r3_melting_conc;
        r3_liquid_new = r3_liquid - r3_melting_conc;
        r3_solid_new = 0.0;
      }
      
    }    
  }
  else
  {
    r3_solid_new = r3_solid;
    r3_liquid_new = r3_liquid;
    r3_gas_new = r3_gas;
  }

  _react1_solid[_qp] = r1_solid_new;
  _react1_liquid[_qp] = r1_liquid_new;
  _react1_gas[_qp] = r1_gas_new;

  _react2_solid[_qp] = r2_solid_new;
  _react2_liquid[_qp] = r2_liquid_new;
  _react2_gas[_qp] = r2_gas_new;

  _react3_solid[_qp] = r3_solid_new;
  _react3_liquid[_qp] = r3_liquid_new;
  _react3_gas[_qp] = r3_gas_new;

  _react4_solid[_qp] = r4_solid_new;
  _react4_liquid[_qp] = r4_liquid_new;
  _react4_gas[_qp] = r4_gas_new;

  // Update the used heat term
  _heat_sink_melting[_qp] = r1_melting_eff + r2_melting_eff + r3_melting_eff + r4_melting_eff;
  _heat_sink_boiling[_qp] = r1_boiling_eff + r2_boiling_eff + r3_boiling_eff + r4_boiling_eff + r4_boiling_eff2;


  // Reactant Phase
  // solid -> 1; solid + liquid -> 1.5; liquid -> 2; liquid + gas -> 2.5; gas -> 3; total reacted -> -1;
  // Reactant 1
  if (r1_solid_new > 0 && r1_liquid_new == 0 && r1_gas_new == 0)
    _react1_phase[_qp] = 1;
  else if (r1_solid_new > 0 && r1_liquid_new > 0 && r1_gas_new == 0)
    _react1_phase[_qp] = 1.5;
  else if (r1_solid_new == 0 && r1_liquid_new > 0 && r1_gas_new == 0)
    _react1_phase[_qp] = 2;
  else if (r1_solid_new ==  0 && r1_liquid_new > 0 && r1_gas_new > 0)
    _react1_phase[_qp] = 2.5;
  else if (r1_solid_new == 0 && r1_liquid_new == 0 && r1_gas_new > 0)
    _react1_phase[_qp] = 3;
  else if (r1_solid_new == 0 && r1_liquid_new == 0 && r1_gas_new == 0)
    _react1_phase[_qp] = 0;
  else
    _react1_phase[_qp] = 3;

  // Reactant 2

  if (r2_solid_new > 0 && r2_liquid_new == 0 && r2_gas_new == 0)
    _react2_phase[_qp] = 1;
  else if (r2_solid_new > 0 && r2_liquid_new > 0 && r2_gas_new == 0)
    _react2_phase[_qp] = 1.5;
  else if (r2_solid_new == 0 && r2_liquid_new > 0 && r2_gas_new == 0)
    _react2_phase[_qp] = 2;
  else if (r2_solid_new ==  0 && r2_liquid_new > 0 && r2_gas_new > 0)
    _react2_phase[_qp] = 2.5;
  else if (r2_solid_new == 0 && r2_liquid_new == 0 && r2_gas_new > 0)
    _react2_phase[_qp] = 3;
  else if (r2_solid_new == 0 && r2_liquid_new == 0 && r2_gas_new == 0)
    _react2_phase[_qp] = 0;
  else
    _react2_phase[_qp] = 3;


  // Reactant 3
  if (r3_solid_new > 0 && r3_liquid_new == 0 && r3_gas_new == 0)
    _react3_phase[_qp] = 1;
  else if (r3_solid_new > 0 && r3_liquid_new > 0 && r3_gas_new == 0)
    _react3_phase[_qp] = 1.5;
  else if (r3_solid_new == 0 && r3_liquid_new > 0 && r3_gas_new == 0)
    _react3_phase[_qp] = 2;
  else if (r3_solid_new ==  0 && r3_liquid_new > 0 && r3_gas_new > 0)
    _react3_phase[_qp] = 2.5;
  else if (r3_solid_new == 0 && r3_liquid_new == 0 && r3_gas_new > 0)
    _react3_phase[_qp] = 3;
  else if (r3_solid_new == 0 && r3_liquid_new == 0 && r3_gas_new == 0)
    _react3_phase[_qp] = 0;
  else
    _react3_phase[_qp] = 3;


  if (r4_solid_new > 0 && r4_liquid_new == 0 && r4_gas_new == 0)
    _react4_phase[_qp] = 1;
  else if (r4_solid_new > 0 && r4_liquid_new > 0 && r4_gas_new == 0)
    _react4_phase[_qp] = 1.5;
  else if (r4_solid_new == 0 && r4_liquid_new > 0 && r4_gas_new == 0)
    _react4_phase[_qp] = 2;
  else if (r4_solid_new ==  0 && r4_liquid_new > 0 && r4_gas_new > 0)
    _react4_phase[_qp] = 2.5;
  else if (r4_solid_new == 0 && r4_liquid_new == 0 && r4_gas_new > 0)
    _react4_phase[_qp] = 3;
  else if (r4_solid_new == 0 && r4_liquid_new == 0 && r4_gas_new == 0)
    _react4_phase[_qp] = 0;
  else
    _react4_phase[_qp] = 3;

  _react1_status[_qp] = _react1_solid[_qp] + _react1_liquid[_qp] + _react1_gas[_qp];
  _react2_status[_qp] = _react2_solid[_qp] + _react2_liquid[_qp] + _react2_gas[_qp];
  _prod1_status[_qp] = _react3_solid[_qp] + _react3_liquid[_qp] + _react3_gas[_qp];
  _prod2_status[_qp] = _react4_solid[_qp] + _react4_liquid[_qp] + _react4_gas[_qp];
  
// Update all different properties, density, thermal conductivity and Cp
// According to the mass fraction of reactants and products
  Real km1(0.0);
  Real km2(0.0);

  for (unsigned int i = 0; i < _species.size(); ++i)
  {
    km1 += _current_vol_fraction[i]/_component_thermal_cond[i];
    km2 += _current_vol_fraction[i]*_component_thermal_cond[i];  
  }
  
  _thermal_conductivity[_qp] = (1/km1+km2)/2;
  
  
  // for now the thermal conductivity and density are treated as constants
  _density[_qp] = 1/rho;

  // loop over all species, check phases
  std::vector<Real> Cp;
  Cp.resize(_species.size());
    
  for (unsigned int i = 0; i < _species.size(); ++i)
  {
    Cp[i] = _specific_heat_solid[i];
    
    if ( T >= _melting_temp[i] && T < _vaporization_temp[i])
      Cp[i] = _specific_heat_liquid[i];

    if ( T >=  _vaporization_temp[i])
      Cp[i] = _specific_heat_liquid[i];
  }
  
  Real effective_Cp(0.0);
    
  for (unsigned int i = 0; i < _species.size(); ++i)
    effective_Cp += Cp[i] * _current_mass_conc[i] / total_mass_conc;
    
  _specific_heat[_qp] = effective_Cp;
}
