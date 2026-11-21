# two-plastic-phase input (biphasic_benchmark, writer r1, pin_mode = corners). Plane-stress periodic RVE, stress control.
[GlobalParams]
  displacements = 'disp_x disp_y'
  large_kinematics = false
  stabilize_strain = false
  macro_gradient = hvar
  constraint_types = 'stress none none stress stress none none none none'
  targets = 'sxx_t sxy_t syy_t'
[]

[Mesh]
  use_displaced_mesh = false
  [file]
    type = FileMeshGenerator
    file = '/home/jarettpoliner/projects/yield_surface_moose/plane_stress/specimens/S2/mesh_h40/rve.msh'
  []
  [pin]
    type = ExtraNodesetGenerator
    input = file
    coord = '0 0 0; 1.0 0 0; 0 1.0 0; 1.0 1.0 0'
    new_boundary = pin
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [strain_zz]
    family = L2_LAGRANGE
    order = FIRST
  []
  [hvar]
    family = SCALAR
    order = THIRD
  []
[]

[Functions]
  [sxx_t]
    type = ParsedFunction
    expression = '1*t'
  []
  [sxy_t]
    type = ParsedFunction
    expression = '0*t'
  []
  [syy_t]
    type = ParsedFunction
    expression = '0*t'
  []
[]

[Kernels]
  [sdx]
    type = HomogenizedTotalLagrangianStressDivergence
    variable = disp_x
    component = 0
    macro_var = hvar
    out_of_plane_strain = strain_zz
  []
  [sdy]
    type = HomogenizedTotalLagrangianStressDivergence
    variable = disp_y
    component = 1
    macro_var = hvar
    out_of_plane_strain = strain_zz
  []
  [wps]
    type = TotalLagrangianWeakPlaneStress
    variable = strain_zz
  []
[]

[BCs]
  [Periodic]
    [x]
      variable = 'disp_x disp_y'
      primary = left
      secondary = right
      translation = '1.0 0 0'
    []
    [y]
      variable = 'disp_x disp_y'
      primary = bottom
      secondary = top
      translation = '0 1.0 0'
    []
  []
  [pinx]
    type = DirichletBC
    boundary = pin
    variable = disp_x
    value = 0
  []
  [piny]
    type = DirichletBC
    boundary = pin
    variable = disp_y
    value = 0
  []
[]

[Materials]
  [strain]
    type = ComputeLagrangianWPSStrain
    out_of_plane_strain = strain_zz
    homogenization_gradient_names = homogenization_gradient
  []
  [hgrad]
    type = ComputeHomogenizedLagrangianStrain
  []
  [C_matrix]
    type = ComputeIsotropicElasticityTensor
    block = matrix
    youngs_modulus = 72000.0
    poissons_ratio = 0.34
  []
  [stress_matrix]
    type = ComputeLagrangianWrappedStress
    block = matrix
    rotate_old_stress = false
  []
  [stress_matrix_base]
    type = ComputeMultiPlasticityStress
    block = matrix
    plastic_models = j2_m
    ep_plastic_tolerance = 1e-10
    tangent_operator = nonlinear
  []
  [C_inclusion]
    type = ComputeIsotropicElasticityTensor
    block = inclusion
    youngs_modulus = 193000.0
    poissons_ratio = 0.29
  []
  [stress_inclusion]
    type = ComputeLagrangianWrappedStress
    block = inclusion
    rotate_old_stress = false
  []
  [stress_inclusion_base]
    type = ComputeMultiPlasticityStress
    block = inclusion
    plastic_models = j2_i
    ep_plastic_tolerance = 1e-10
    tangent_operator = nonlinear
  []
  [sxx]
    type = RankTwoCartesianComponent
    rank_two_tensor = cauchy_stress
    property_name = sxx
    index_i = 0
    index_j = 0
  []
  [syy]
    type = RankTwoCartesianComponent
    rank_two_tensor = cauchy_stress
    property_name = syy
    index_i = 1
    index_j = 1
  []
  [sxy]
    type = RankTwoCartesianComponent
    rank_two_tensor = cauchy_stress
    property_name = sxy
    index_i = 0
    index_j = 1
  []
  [szz]
    type = RankTwoCartesianComponent
    rank_two_tensor = cauchy_stress
    property_name = szz
    index_i = 2
    index_j = 2
  []
  [sxz]
    type = RankTwoCartesianComponent
    rank_two_tensor = cauchy_stress
    property_name = sxz
    index_i = 0
    index_j = 2
  []
  [syz]
    type = RankTwoCartesianComponent
    rank_two_tensor = cauchy_stress
    property_name = syz
    index_i = 1
    index_j = 2
  []
  [exx]
    type = RankTwoCartesianComponent
    rank_two_tensor = total_strain
    property_name = exx
    index_i = 0
    index_j = 0
  []
  [eyy]
    type = RankTwoCartesianComponent
    rank_two_tensor = total_strain
    property_name = eyy
    index_i = 1
    index_j = 1
  []
  [exy]
    type = RankTwoCartesianComponent
    rank_two_tensor = total_strain
    property_name = exy
    index_i = 0
    index_j = 1
  []
  [ezz]
    type = RankTwoCartesianComponent
    rank_two_tensor = total_strain
    property_name = ezz
    index_i = 2
    index_j = 2
  []
  [vm]
    type = RankTwoInvariant
    rank_two_tensor = cauchy_stress
    property_name = vm
    invariant = VonMisesStress
  []
[]

[Postprocessors]
  [s_xx]
    type = ElementAverageMaterialProperty
    mat_prop = sxx
    execute_on = 'initial timestep_end'
  []
  [e_xx]
    type = ElementAverageMaterialProperty
    mat_prop = exx
    execute_on = 'initial timestep_end'
  []
  [s_yy]
    type = ElementAverageMaterialProperty
    mat_prop = syy
    execute_on = 'initial timestep_end'
  []
  [e_yy]
    type = ElementAverageMaterialProperty
    mat_prop = eyy
    execute_on = 'initial timestep_end'
  []
  [s_xy]
    type = ElementAverageMaterialProperty
    mat_prop = sxy
    execute_on = 'initial timestep_end'
  []
  [e_xy]
    type = ElementAverageMaterialProperty
    mat_prop = exy
    execute_on = 'initial timestep_end'
  []
  [s_zz]
    type = ElementAverageMaterialProperty
    mat_prop = szz
    execute_on = 'initial timestep_end'
  []
  [e_zz]
    type = ElementAverageMaterialProperty
    mat_prop = ezz
    execute_on = 'initial timestep_end'
  []
  [szz_qp_max]
    type = ElementExtremeMaterialProperty
    mat_prop = szz
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [szz_qp_min]
    type = ElementExtremeMaterialProperty
    mat_prop = szz
    value_type = min
    execute_on = 'initial timestep_end'
  []
  [sxz_qp_max]
    type = ElementExtremeMaterialProperty
    mat_prop = sxz
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [sxz_qp_min]
    type = ElementExtremeMaterialProperty
    mat_prop = sxz
    value_type = min
    execute_on = 'initial timestep_end'
  []
  [syz_qp_max]
    type = ElementExtremeMaterialProperty
    mat_prop = syz
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [syz_qp_min]
    type = ElementExtremeMaterialProperty
    mat_prop = syz
    value_type = min
    execute_on = 'initial timestep_end'
  []
  [eqps_rve_avg]
    type = ElementAverageValue
    variable = eqps
    execute_on = 'initial timestep_end'
  []
  [wp_rve_avg]
    type = ElementAverageValue
    variable = wp
    execute_on = 'initial timestep_end'
  []
  [eqps_m_avg]
    type = ElementAverageValue
    variable = eqps
    block = matrix
    execute_on = 'initial timestep_end'
  []
  [eqps_m_max]
    type = ElementExtremeValue
    variable = eqps
    value_type = max
    block = matrix
    execute_on = 'initial timestep_end'
  []
  [pvf_m]
    type = ElementAverageValue
    variable = plastic_ind
    block = matrix
    execute_on = 'initial timestep_end'
  []
  [wp_m_avg]
    type = ElementAverageValue
    variable = wp
    block = matrix
    execute_on = 'initial timestep_end'
  []
  [vm_m_max]
    type = ElementExtremeMaterialProperty
    mat_prop = vm
    value_type = max
    block = matrix
    execute_on = 'initial timestep_end'
  []
  [vm_m_avg]
    type = ElementAverageMaterialProperty
    mat_prop = vm
    block = matrix
    execute_on = 'initial timestep_end'
  []
  [s_xx_m_avg]
    type = ElementAverageMaterialProperty
    mat_prop = sxx
    block = matrix
    execute_on = 'initial timestep_end'
  []
  [s_yy_m_avg]
    type = ElementAverageMaterialProperty
    mat_prop = syy
    block = matrix
    execute_on = 'initial timestep_end'
  []
  [s_xy_m_avg]
    type = ElementAverageMaterialProperty
    mat_prop = sxy
    block = matrix
    execute_on = 'initial timestep_end'
  []
  [vol_m]
    type = VolumePostprocessor
    block = matrix
    execute_on = initial
  []
  [eqps_i_avg]
    type = ElementAverageValue
    variable = eqps
    block = inclusion
    execute_on = 'initial timestep_end'
  []
  [eqps_i_max]
    type = ElementExtremeValue
    variable = eqps
    value_type = max
    block = inclusion
    execute_on = 'initial timestep_end'
  []
  [pvf_i]
    type = ElementAverageValue
    variable = plastic_ind
    block = inclusion
    execute_on = 'initial timestep_end'
  []
  [wp_i_avg]
    type = ElementAverageValue
    variable = wp
    block = inclusion
    execute_on = 'initial timestep_end'
  []
  [vm_i_max]
    type = ElementExtremeMaterialProperty
    mat_prop = vm
    value_type = max
    block = inclusion
    execute_on = 'initial timestep_end'
  []
  [vm_i_avg]
    type = ElementAverageMaterialProperty
    mat_prop = vm
    block = inclusion
    execute_on = 'initial timestep_end'
  []
  [s_xx_i_avg]
    type = ElementAverageMaterialProperty
    mat_prop = sxx
    block = inclusion
    execute_on = 'initial timestep_end'
  []
  [s_yy_i_avg]
    type = ElementAverageMaterialProperty
    mat_prop = syy
    block = inclusion
    execute_on = 'initial timestep_end'
  []
  [s_xy_i_avg]
    type = ElementAverageMaterialProperty
    mat_prop = sxy
    block = inclusion
    execute_on = 'initial timestep_end'
  []
  [vol_i]
    type = VolumePostprocessor
    block = inclusion
    execute_on = initial
  []
  [vm_max]
    type = ElementExtremeMaterialProperty
    mat_prop = vm
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [nl_its]
    type = NumNonlinearIterations
    execute_on = timestep_end
  []
  [l_its]
    type = NumLinearIterations
    execute_on = timestep_end
  []
[]

[AuxVariables]
  [eqps]
    family = MONOMIAL
    order = CONSTANT
  []
  [plastic_ind]
    family = MONOMIAL
    order = CONSTANT
  []
  [wp]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [eqps]
    type = MaterialStdVectorAux
    variable = eqps
    property = plastic_internal_parameter
    index = 0
    block = 'matrix inclusion'
    execute_on = 'initial timestep_end'
  []
  [plastic_ind]
    type = ParsedAux
    variable = plastic_ind
    coupled_variables = eqps
    expression = 'if(eqps > 1e-06, 1, 0)'
    execute_on = 'initial timestep_end'
  []
  [wp_m]
    type = ParsedAux
    variable = wp
    block = matrix
    coupled_variables = eqps
    expression = '200*eqps + 0.5*800*eqps*eqps'
    execute_on = 'initial timestep_end'
  []
  [wp_i]
    type = ParsedAux
    variable = wp
    block = inclusion
    coupled_variables = eqps
    expression = '215*eqps + 0.5*1500*eqps*eqps'
    execute_on = 'initial timestep_end'
  []
[]

[UserObjects]
  [hardening_m]
    type = SolidMechanicsHardeningPowerRule
    value_0 = 200
    epsilon0 = 0.25
    exponent = 1
  []
  [j2_m]
    type = SolidMechanicsPlasticJ2
    yield_strength = hardening_m
    yield_function_tolerance = 1e-6
    internal_constraint_tolerance = 1e-10
  []
  [hardening_i]
    type = SolidMechanicsHardeningPowerRule
    value_0 = 215
    epsilon0 = 0.1433333333
    exponent = 1
  []
  [j2_i]
    type = SolidMechanicsPlasticJ2
    yield_strength = hardening_i
    yield_function_tolerance = 1e-6
    internal_constraint_tolerance = 1e-10
  []
  [stop]
    type = Terminator
    expression = 'eqps_rve_avg > 0.004'
    fail_mode = HARD
    execute_on = timestep_end
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu superlu_dist'
  line_search = none
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-9
  nl_max_its = 30
  l_tol = 1e-8
  l_max_its = 50
  [TimeStepper]
    type = TimeSequenceStepper
    time_sequence = '1 137.2801308 144.1441373 151.3513442 158.9189114 166.8648569 175.2080998 183.9685048 193.16693 202.8252765 212.9665403 223.6148674 234.7956107 246.5353913 258.8621608 271.8052689 285.3955323 299.6653089 314.6485744 330.3810031 346.9000533 364.2450559 382.4573087 401.5801741 421.6591829 442.742142 464.8792491 488.1232116 512.5293721 538.1558407 565.0636328 593.3168144 622.9826551 654.1317879 686.8383773 721.1802962 757.239311 795.1012765 834.8563403 876.5991573 920.4291152'
  []
  start_time = 0
  end_time = 920.4291152
[]

[Outputs]
  print_linear_residuals = false
  [csv]
    type = CSV
    file_base = 'ramp_out'
    execute_vector_postprocessors_on = final
  []
  [cp]
    type = Checkpoint
    num_files = 52
  []
[]
