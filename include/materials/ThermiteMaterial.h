#ifndef THERMITEMATERIAL_H
#define THERMITEMATERIAL_H

#include "Material.h"

class ThermiteMaterial;

template<>
InputParameters validParams<ThermiteMaterial>();

class ThermiteMaterial : public Material
{
public:
  ThermiteMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

private:

  Real _initial_thermal_conductivity;
  Real _initial_specific_heat;
  Real _initial_density;

  Real _ignition_temp;
  Real _kr;
  

  MaterialProperty<Real> & _thermal_conductivity; //thermal_conductivity
  MaterialProperty<Real> & _specific_heat; //specific heat
  MaterialProperty<Real> & _density; //mixutre density
  
  MaterialProperty<Real> & _r_ex;     //reaction extent
  MaterialProperty<Real> & _r_ex_old; //reaction extent at previous time step
  
  const VariableValue & _temp; //temp

  ///thermodynamic properties of thermite
  std::vector<std::string> _species;
  std::vector<Real> _sto_coeff;
  std::vector<Real> _molecular_weights;
  std::vector<Real> _melting_temp;
  std::vector<Real> _vaporization_temp;
  std::vector<Real> _latent_s_to_l;
  std::vector<Real> _latent_l_to_g;
  

  std::vector<Real> _initial_mass_fraction;
  std::vector<Real> _current_mass_fraction;
  std::vector<Real> _initial_mass_conc; // kg/m^3
  std::vector<Real> _current_mass_conc; // kg/m^3


  std::vector<Real> _initial_molar_conc; // mols/m^3
  std::vector<Real> _current_molar_conc; // mols/m^3
  

  
  

};

#endif //THERMITEMATERIAL_H
