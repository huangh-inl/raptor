/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef COMBUST_H
#define COMBUST_H

#include "Kernel.h"
#include "Reacted.h"

//Forward Declarations
class combust;
class Function;

template<>
InputParameters validParams<combust>();

/**
 * Define the Kernel for a user defined forcing function that looks like:
 *
 * test function * forcing function
 */
class combust : public Kernel
{
public:

  combust(const InputParameters & parameters);

protected:
  /**
   * Evaluate f at the current quadrature point.
   */
  Real f();

  /**
   * Computes test function * forcing function.
   */
  virtual Real computeQpResidual();

  const MaterialProperty<Real> & _r_ex;
  const MaterialProperty<Real> & _r_ex_old;
  const MaterialProperty<Real> & _heat_sink_melting;
  const MaterialProperty<Real> & _heat_sink_boiling;
  
//   const MaterialProperty<Real> & _thermal_conductivity;
//   const MaterialProperty<Real> & _specific_heat;
//   const MaterialProperty<Real> & _density;
//   const VariableValue & _coupled_val;
//   const VariableValue & _dcoupled_val_dt;

  Function & _func;
};

#endif //COMBUST_H
