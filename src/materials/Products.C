#include "Products.h"

template<>
InputParameters validParams<Products>()
{
  InputParameters params = validParams<Material>();

  params.addCoupledVar("coupled", "Status of the reaction");
  params.addCoupledVar("temperature", "Temperature at current point");
  params.addParam<Real>("initial_thermal_conductivity", 1.0, "The diffusion coefficient of the Reactants");
  params.addParam<Real>("initial_specific_heat", 1.0, "The specific heat of the Reactants");
  params.addParam<Real>("initial_density", 1.0, "The density the Reactants");

  return params;
}

Products::Products(const InputParameters & parameters) :
  Material(parameters),

  // Get a parameter value for the specific heat
  _initial_thermal_conductivity(getParam<Real>("initial_thermal_conductivity")),
  _initial_specific_heat(getParam<Real>("initial_specific_heat")),
  _initial_density(getParam<Real>("initial_density")),

  // Declare that this material is going to have a real
  // value property named "specific_heat" that Kernels can use.
  _thermal_conductivity(declareProperty<Real>("thermal_conductivity")),  
  _specific_heat(declareProperty<Real>("specific_heat")),  
  _density(declareProperty<Real>("density")),  

  _coupled_status(coupledValue("coupled")),
  _temp(coupledValue("temperature")) 
{}

// void
// Products::initQpStatefulProperties()
// {
  // init the diffusivity property (this will become
  // _diffusivity_old in the first call of computeProperties)
//   _specific_heat[_qp] = _initial_specific_heat;
// }


void
Products::computeQpProperties()
{
  if (_coupled_status[_qp] > 0.0 && _temp[_qp] > 1100)
  {
    _thermal_conductivity[_qp] =  _initial_thermal_conductivity * 0.5;
    _specific_heat[_qp] = _initial_specific_heat / 0.89;
    _density[_qp] = _initial_density * 0.768;
  }
  else
  {
    _thermal_conductivity[_qp] =  _initial_thermal_conductivity;
    _specific_heat[_qp] = _initial_specific_heat;
    _density[_qp] = _initial_density;
    
//    _thermal_conductivity[_qp] = _temp[_qp] / 100.0;
    //_thermal_conductivity[_qp] = 20;
//    _specific_heat[_qp] = 880;
//    _density[_qp] = 4500;
  }
//  std::cout<<_thermal_conductivity[_qp]<<std::endl;
}
