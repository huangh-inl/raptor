[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 400
  ny = 1
#  xmax = 0.0258 # Length of test chamber
  xmax = 0.1032
#   xmax = 0.00387
  ymax = 0.000258 # Test chamber radius
[]

[Variables]
  [./temp]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = BoundingBoxIC
      x1 = 0.0
      y1 = 0.0
      x2 = 0.00258
      y2 = 0.0108
      inside = 1210
      outside = 300
    [../]

#    initial_condition = 300 # Start at room temperature
  [../]
[]

[AuxVariables]
  [./reaction_extent]
    order = CONSTANT
    family = Monomial
  [../]
  # [./thermal_conductivity]
  #   order = CONSTANT
  #   family = Monomial
  # [../]
  # [./specific_heat]
  #   order = CONSTANT
  #   family = Monomial
  # [../]
  [./Al_phase]
    order = CONSTANT
    family = Monomial
  [../]
  [./Fe2O3_phase]
    order = CONSTANT
    family = Monomial
  [../]
  [./Al2O3_phase]
    order = CONSTANT
    family = Monomial
  [../]
  [./Fe_phase]
    order = CONSTANT
    family = Monomial
  [../]
  [./Al_state]
    order = CONSTANT
    family = Monomial
  [../]
  [./Fe2O3_state]
    order = CONSTANT
    family = Monomial
  [../]
  [./Al2O3_state]
    order = CONSTANT
    family = Monomial
  [../]
  [./Fe_state]
    order = CONSTANT
    family = Monomial
  [../]
  [./Al_solid]
    order = CONSTANT
    family = Monomial
  [../]
  [./Al_liquid]
    order = CONSTANT
    family = Monomial
  [../]
  [./Al_gas]
    order = CONSTANT
    family = Monomial
  [../]
  [./Fe2O3_solid]
    order = CONSTANT
    family = Monomial
  [../]
  [./Fe2O3_liquid]
    order = CONSTANT
    family = Monomial
  [../]
  [./Fe2O3_gas]
    order = CONSTANT
    family = Monomial
  [../]
  [./Al2O3_solid]
    order = CONSTANT
    family = Monomial
  [../]
  [./Al2O3_liquid]
    order = CONSTANT
    family = Monomial
  [../]
  [./Al2O3_gas]
    order = CONSTANT
    family = Monomial
  [../]
  [./Fe_solid]
    order = CONSTANT
    family = Monomial
  [../]
  [./Fe_liquid]
    order = CONSTANT
    family = Monomial
  [../]
  [./Fe_gas]
    order = CONSTANT
    family = Monomial
  [../]
        
 []

[Functions]
  [./Combust]
    type = ParsedFunction
    value = 5322.746*alpha
    vars = alpha
    vals =  4360
  [../]
[]

[Kernels]
  [./heat_conduction]
    type = HeatConduction
    variable = temp
  [../]
  [./heat_conduction_time_derivative]
    type = HeatConductionTimeDerivative
    variable = temp
  [../]

  [./reaction_heat_source]
   type = combust
   variable = temp
   function = Combust
#   coupled = Reaction_extend
  [../]
[]

[AuxKernels]
  [./reaction_status]
    type = MaterialRealAux
    variable = reaction_extent
    property = reaction_extent
    execute_on = 'initial timestep_end'
  [../]
  # [./thermal_k]
  #   type = MaterialRealAux
  #   variable = thermal_conductivity
  #   property = thermal_conductivity
  #   execute_on = 'initial timestep_end'
  # [../]
  # [./cp]
  #   type = MaterialRealAux
  #   variable = specific_heat
  #   property = specific_heat
  #   execute_on = 'initial timestep_end'
  # [../]
  [./react1]
    type = MaterialRealAux
    variable = Al_phase
    property = Al_phase
    execute_on = 'initial timestep_end'
  [../]
  [./react2]
    type = MaterialRealAux
    variable = Fe2O3_phase
    property = Fe2O3_phase
    execute_on = 'initial timestep_end'
  [../]
  [./product1]
    type = MaterialRealAux
    variable = Al2O3_phase
    property = Al2O3_phase
    execute_on = 'initial timestep_end'
  [../]
  [./product2]
    type = MaterialRealAux
    variable = Fe_phase
    property = Fe_phase
    execute_on = 'initial timestep_end'
  [../]
  [./react11]
    type = MaterialRealAux
    variable = Al_state
    property = Al_state
    execute_on = 'initial timestep_end'
  [../]
  [./react21]
    type = MaterialRealAux
    variable = Fe2O3_state
    property = Fe2O3_state
    execute_on = 'initial timestep_end'
  [../]
  [./product11]
    type = MaterialRealAux
    variable = Al2O3_state
    property = Al2O3_state
    execute_on = 'initial timestep_end'
  [../]
  [./product21]
    type = MaterialRealAux
    variable = Fe_state
    property = Fe_state
    execute_on = 'initial timestep_end'
  [../]
  [./Al1]
    type = MaterialRealAux
    variable = Al_solid
    property = Al_solid
    execute_on = 'initial timestep_end'
  [../]
  [./Al2]
    type = MaterialRealAux
    variable = Al_liquid
    property = Al_liquid
    execute_on = 'initial timestep_end'
  [../]
  [./Al3]
    type = MaterialRealAux
    variable = Al_gas
    property = Al_gas
    execute_on = 'initial timestep_end'
  [../]
  [./Fe2O31]
    type = MaterialRealAux
    variable = Fe2O3_solid
    property = Fe2O3_solid
    execute_on = 'initial timestep_end'
  [../]
  [./Fe2O32]
    type = MaterialRealAux
    variable = Fe2O3_liquid
    property = Fe2O3_liquid
    execute_on = 'initial timestep_end'
  [../]
  [./Fe2O33]
    type = MaterialRealAux
    variable = Fe2O3_gas
    property = Fe2O3_gas
    execute_on = 'initial timestep_end'
  [../]
  [./Al2O31]
    type = MaterialRealAux
    variable = Al2O3_solid
    property = Al2O3_solid
    execute_on = 'initial timestep_end'
  [../]
  [./Al2O32]
    type = MaterialRealAux
    variable = Al2O3_liquid
    property = Al2O3_liquid
    execute_on = 'initial timestep_end'
  [../]
  [./Al2O33]
    type = MaterialRealAux
    variable = Al2O3_gas
    property = Al2O3_gas
    execute_on = 'initial timestep_end'
  [../]
  [./Fe1]
    type = MaterialRealAux
    variable = Fe_solid
    property = Fe_solid
    execute_on = 'initial timestep_end'
  [../]
  [./Fe2]
    type = MaterialRealAux
    variable = Fe_liquid
    property = Fe_liquid
    execute_on = 'initial timestep_end'
  [../]
  [./Fe3]
    type = MaterialRealAux
    variable = Fe_gas
    property = Fe_gas
    execute_on = 'initial timestep_end'
  [../]
[]


[BCs]
  #[./input_temperature]
  #  type = DirichletBC
  #  variable = temp
  #  boundary = left
  #  value = 2390
  #[../]
  [./inlet_temperature]
    type = ConvectiveFluxBC
    variable = temp
    final = 300.0
    rate = 1.0
    boundary = left
  [../]
  [./outlet_temperature]
    type = ConvectiveFluxBC
    variable = temp
    final = 300.0
    rate = 1.0
    boundary = right
  [../]
  #[./output_temperature]
  #  type = NeumannBC
  #  variable = temp
  #  boundary = right
  #  function = Combust
  #[../]
[]

[Materials]
  active = 'reactants'
  [./reactants]
    type = ThermiteMaterial
    initial_thermal_conductivity = 10.0e-3  # kW/m*K
    initial_specific_heat = 0.680  # kJ/kg-K
    initial_density = 4360  # kg/m^3
    heat_fraction = 0.5

    temp_var = temp
    reaction_system =   'Al Fe2O3 Al2O3 Fe'
    sto_coeff =         '2   1      1   2'
    molecular_weights = '26.98e-3 159.69e-3   101.96e-3  55.85e-3'
    melting_temp =      '933.15   1839.15     2369.0     1808.0'
    vaporization_temp = '2740.15  3250.0      5950.15    3135.15'
    latent_s_to_l     = '398.00   126.0       1360.0     272.0'
    latent_l_to_g     = '11400.0  3088.0      35000.0    6087.74'
    specific_heat_solid = '0.91  0.78   0.955  0.45'
    specific_heat_liquid = '1.18 1.35   1.88   0.82'
    component_density = '2700.0  5240.0 4000.0 7870.0'
    component_thermal_cond = '0.205  0.01  0.0385  0.0795'

    ignition_temp     = 1200.0
    rate_constant = 1.0
    
  [../]
[]

# [Problem]
#   type = FEProblem
#   coord_type = RZ
#   rz_coord_axis = X
# []

[Executioner]
  type       = Transient
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -snes_ls -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 201 cubic 0.7'


  dtmax = 0.1
  dtmin = 1.0e-12
  end_time = 0.001 # 0.001 for test | 11 for full run

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 1.0e-5
  [../]

  l_max_its  = 50
  l_tol      = 1e-5
  nl_max_its = 20
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
[]

[Outputs]
  file_base      = heat05
  exodus         = true
  interval       = 10
  [./Console]
    type = Console
    linear_residuals    = 1
    nonlinear_residuals = 1
  [../]
[]
