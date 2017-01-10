#include "combust.h"
#include "Function.h"

template<>
InputParameters validParams<combust>()
{
  InputParameters params = validParams<Kernel>();
//  params.addRequiredCoupledVar("coupled", "Coupled variable (reaction status)");
  params.addRequiredParam<FunctionName>("function", "The forcing function");
  return params;
}

combust::combust(const InputParameters & parameters) :
    Kernel(parameters),
//     _thermal_conductivity(getMaterialProperty<Real>("thermal_conductivity")),
//     _specific_heat(getMaterialProperty<Real>("specific_heat")),
//     _density(getMaterialProperty<Real>("density")),
//     _coupled_val(coupledValue("coupled")),
//     _dcoupled_val_dt(coupledDot("coupled")),
    _r_ex(getMaterialProperty<Real>("reaction_extent")),
    _r_ex_old(getMaterialPropertyOld<Real>("reaction_extent")),
    _heat_sink_melting(getMaterialProperty<Real>("melting_heatsink")),
    _heat_sink_boiling(getMaterialProperty<Real>("boiling_heatsink")),
    _func(getFunction("function"))
{}

Real
combust::f()
{
  return _func.value(_t, _q_point[_qp]);
}

Real
combust::computeQpResidual()
{

//   if(_dcoupled_val_dt[_qp] > 0.0 && _coupled_val[_qp] <= 1.0)
    
//     return -_test[_i][_qp] * f() * _dcoupled_val_dt[_qp];
//   else
//     return 0.0;

  // f() is the total reaction heat KJ/m^3 = ro/ave_molecular_weight * dH, where dH is reaction heat in KJ/mole
  Real dval_dt = (_r_ex[_qp] - _r_ex_old[_qp])/_dt;
  Real phase_heat = _heat_sink_melting[_qp] + _heat_sink_boiling[_qp];
  Real tot_heat = f() * dval_dt - phase_heat;
  Real tot_r_heat = f() * dval_dt;
  

  if (tot_heat < 0)
    tot_heat = 0;

  // if (tot_r_heat > 0)
  //   std::cout << tot_r_heat << "  " << phase_heat << "  "  << tot_heat << std::endl;
  
  return -_test[_i][_qp] * tot_heat;
  
}
