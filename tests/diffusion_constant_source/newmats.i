[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 1
  xmax = 0.0258 # Length of test chamber
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
  [./thermal_conductivity]
    order = CONSTANT
    family = Monomial
  [../]
  [./specific_heat]
    order = CONSTANT
    family = Monomial
  [../]
[]

[Functions]
  [./Combust]
    type = ParsedFunction
    value = 532200.746*alpha
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
  [./thermal_k]
    type = MaterialRealAux
    variable = thermal_conductivity
    property = thermal_conductivity
    execute_on = 'initial timestep_end'
  [../]
  [./cp]
    type = MaterialRealAux
    variable = specific_heat
    property = specific_heat
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
    initial_thermal_conductivity = 7.0  # kW/m*K
    initial_specific_heat = 0.680  # kJ/kg-K
    initial_density = 4360  # kg/m^3

    temp_var = temp
    reaction_system =   'Al Fe2O3 Al2O3 Fe'
    sto_coeff =         '2   1      1   2'
    molecular_weights = '26.98e-3 159.69e-3  101.96e-3 55.85e-3'
    melting_temp =      '933.45e3   1838.15e3     2345.14  1811.15'
    vaporization_temp = '2743.15  3500.0      3250.15  3135.15'
    latent_s_to_l     = '321.00   100.0        100.0     247.1'
    latent_l_to_g     = '10500.0  12000        12000    6213.07'

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
  end_time = 5.0

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 1.0e-5
  [../]

  l_max_its  = 50
  l_tol      = 1e-5
  nl_max_its = 20
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-12
[]

# [Postprocessors]
#   [./dofs]
#     type = NumDOFs
#   [../]
#   [./integral]
#     type = ElementL2Error
#     variable = temp
#     function = Combust
#   [../]
# []

[Outputs]
  file_base      = out
  exodus         = true
  interval       = 1
  [./Console]
    type = Console
#    output_linear = 1
#    output_nonlinear = 1
    linear_residuals    = 1
    nonlinear_residuals = 1
  [../]
[]
