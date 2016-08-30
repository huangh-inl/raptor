[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 10
  xmax = 0.0258 # Length of test chamber
  ymax = 0.00258 # Test chamber radius
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
      inside = 1200
      outside = 300
    [../]

#    initial_condition = 300 # Start at room temperature
  [../]
[]

[AuxVariables]
  [./Reaction_extend]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]
[]

[Functions]
  [./Combust]
    type = ParsedFunction
    value = 6000.0*alpha
    vars = alpha
    vals =  1.0e6
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
   coupled = Reaction_extend
  [../]
[]

[AuxKernels]
  active = 'reaction_status'
  [./reaction_status]
    type = Reacted
    variable = Reaction_extend
    ignition_temp = 1100
    rate_constant = 1
    temp_var = temp
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
    rate = 10.0
    boundary = left
  [../]
  [./outlet_temperature]
    type = ConvectiveFluxBC
    variable = temp
    final = 300.0
    rate = 10.0
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
    type = Products
    initial_thermal_conductivity = 0.081  # W/m*K
    initial_specific_heat = 0.081  # J/kg-K
    initial_density = 4360  # kg/m^3
    coupled = Reaction_extend
    temperature = temp
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
  dtmin = 1.0e-9
  end_time = 5.0

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 1.0e-5
  [../]

  l_max_its  = 50
  l_tol      = 1e-6
  nl_max_its = 10
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-15
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
  execute_on     = 'timestep_end'
  exodus         = true
  interval       = 1
  [./Console]
    type = Console
    output_linear = 1
    output_nonlinear = 1
    #linear_residuals    = 1
    #nonlinear_residuals = 1
  [../]
[]
