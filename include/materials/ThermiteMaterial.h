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
  Real _heat_fraction;

  Real _ignition_temp;
  Real _kr;
  

  MaterialProperty<Real> & _thermal_conductivity; //thermal_conductivity
  MaterialProperty<Real> & _specific_heat; //specific heat
  MaterialProperty<Real> & _density; //mixutre density
  
  MaterialProperty<Real> & _r_ex;     //reaction extent
  MaterialProperty<Real> & _r_ex_old; //reaction extent at previous time step

  MaterialProperty<Real> & _heat_sink_melting;  // heat sink term due to phase change from solid to liquid
  MaterialProperty<Real> & _heat_sink_boiling;  // heat sink term due to phase change from liquid to vapor

  MaterialProperty<Real> & _react1_status;
  MaterialProperty<Real> & _react2_status;
  MaterialProperty<Real> & _prod1_status;
  MaterialProperty<Real> & _prod2_status;

  MaterialProperty<Real> & _react1_phase;
  MaterialProperty<Real> & _react2_phase;
  MaterialProperty<Real> & _react3_phase;
  MaterialProperty<Real> & _react4_phase;

  MaterialProperty<Real> & _react1_solid;
  MaterialProperty<Real> & _react1_solid_old;
  MaterialProperty<Real> & _react1_liquid;
  MaterialProperty<Real> & _react1_liquid_old;
  MaterialProperty<Real> & _react1_gas;
  MaterialProperty<Real> & _react1_gas_old;

  MaterialProperty<Real> & _react2_solid;
  MaterialProperty<Real> & _react2_solid_old;
  MaterialProperty<Real> & _react2_liquid;
  MaterialProperty<Real> & _react2_liquid_old;
  MaterialProperty<Real> & _react2_gas;
  MaterialProperty<Real> & _react2_gas_old;

  MaterialProperty<Real> & _react3_solid;
  MaterialProperty<Real> & _react3_solid_old;
  MaterialProperty<Real> & _react3_liquid;
  MaterialProperty<Real> & _react3_liquid_old;
  MaterialProperty<Real> & _react3_gas;
  MaterialProperty<Real> & _react3_gas_old;

  MaterialProperty<Real> & _react4_solid;
  MaterialProperty<Real> & _react4_solid_old;
  MaterialProperty<Real> & _react4_liquid;
  MaterialProperty<Real> & _react4_liquid_old;
  MaterialProperty<Real> & _react4_gas;
  MaterialProperty<Real> & _react4_gas_old;
  
  const VariableValue & _temp; //temp

  ///thermodynamic properties of thermite
  std::vector<std::string> _species;
  std::vector<Real> _sto_coeff;
  std::vector<Real> _molecular_weights;
  std::vector<Real> _melting_temp;
  std::vector<Real> _vaporization_temp;
  std::vector<Real> _latent_s_to_l;
  std::vector<Real> _latent_l_to_g;
  std::vector<Real> _specific_heat_solid;
  std::vector<Real> _specific_heat_liquid;
  std::vector<Real> _component_density;
  std::vector<Real> _component_thermal_cond;
  
  std::vector<Real> _initial_mass_fraction;
  std::vector<Real> _current_mass_fraction;
  std::vector<Real> _initial_mass_conc; // kg/m^3
  std::vector<Real> _current_mass_conc; // kg/m^3
  std::vector<Real> _current_vol_fraction;

  std::vector<Real> _initial_molar_conc; // mols/m^3
  std::vector<Real> _current_molar_conc; // mols/m^3

  // Reacted mass concentration of each phase
  std::vector<Real> _reacted_mass_conc; // kg/m^3

};

#endif //THERMITEMATERIAL_H
