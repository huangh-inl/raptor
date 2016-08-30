#ifndef PRODUCTS_H
#define PRODUCTS_H

#include "Material.h"

class Products;

template<>
InputParameters validParams<Products>();

class Products : public Material
{
public:
  Products(const InputParameters & parameters);

protected:
//  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

private:

  Real _initial_thermal_conductivity;
  Real _initial_specific_heat;
  Real _initial_density;
  MaterialProperty<Real> & _thermal_conductivity; //thermal_conductivity
  MaterialProperty<Real> & _specific_heat; //specific heat
  MaterialProperty<Real> & _density; //density
  const VariableValue & _coupled_status; //dx
  const VariableValue & _temp; //temp

};

#endif //PRODUCTS_H
