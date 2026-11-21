# Single-element stress-controlled proportional loading to first yield -- template for MooseYieldOracle.
# UNTESTED against a real MOOSE installation: verify object names for your revision (`<app> --registry`).
# Placeholders filled by the adapter: ${sxx} ${syy} ${sxy} (unit Cauchy direction), ${t_levels} ${t_end} ${run_id},
# and the hidden material parameters ${E} ${nu} ${sigma_y} ${H} (J2 plasticity with linear hardening as the first validation model).
# Time t is the load factor lambda = Mandel radius of the applied stress (MPa); tractions t = lambda * sigma_hat . n on the lateral faces,
# +/-z faces free (plane stress), three-node rigid-body constraints.
[Mesh]
  [cube]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
    elem_type = HEX8
  []
  [n0]
    type = ExtraNodesetGenerator
    input = cube
    new_boundary = pin_xyz
    coord = '0 0 0'
  []
  [n1]
    type = ExtraNodesetGenerator
    input = n0
    new_boundary = pin_yz
    coord = '1 0 0'
  []
  [n2]
    type = ExtraNodesetGenerator
    input = n1
    new_boundary = pin_z
    coord = '0 1 0'
  []
[]
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]
[Physics/SolidMechanics/QuasiStatic/all]
  strain = SMALL
  incremental = true
  add_variables = true
  generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz'
[]
[Functions]
  [txx]
    type = ParsedFunction
    expression = '${sxx}*t'
  []
  [tyy]
    type = ParsedFunction
    expression = '${syy}*t'
  []
  [txy]
    type = ParsedFunction
    expression = '${sxy}*t'
  []
  [mtxx]
    type = ParsedFunction
    expression = '-${sxx}*t'
  []
  [mtyy]
    type = ParsedFunction
    expression = '-${syy}*t'
  []
  [mtxy]
    type = ParsedFunction
    expression = '-${sxy}*t'
  []
[]
[BCs]
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = pin_xyz
    value = 0
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'pin_xyz pin_yz'
    value = 0
  []
  [fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'pin_xyz pin_yz pin_z'
    value = 0
  []
  [right_x]
    type = FunctionNeumannBC
    variable = disp_x
    boundary = right
    function = txx
  []
  [right_y]
    type = FunctionNeumannBC
    variable = disp_y
    boundary = right
    function = txy
  []
  [left_x]
    type = FunctionNeumannBC
    variable = disp_x
    boundary = left
    function = mtxx
  []
  [left_y]
    type = FunctionNeumannBC
    variable = disp_y
    boundary = left
    function = mtxy
  []
  [top_x]
    type = FunctionNeumannBC
    variable = disp_x
    boundary = top
    function = txy
  []
  [top_y]
    type = FunctionNeumannBC
    variable = disp_y
    boundary = top
    function = tyy
  []
  [bottom_x]
    type = FunctionNeumannBC
    variable = disp_x
    boundary = bottom
    function = mtxy
  []
  [bottom_y]
    type = FunctionNeumannBC
    variable = disp_y
    boundary = bottom
    function = mtyy
  []
[]
[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = j2
  []
  [j2]
    type = IsotropicPlasticityStressUpdate
    yield_stress = ${sigma_y}
    hardening_constant = ${H}
  []
[]
[AuxVariables]
  [eqps_aux]
    order = CONSTANT
    family = MONOMIAL
  []
[]
[AuxKernels]
  [eqps_aux]
    type = MaterialRealAux
    property = effective_plastic_strain
    variable = eqps_aux
  []
[]
[Postprocessors]
  [eqps]
    type = ElementAverageValue
    variable = eqps_aux
  []
  [sxx]
    type = ElementAverageValue
    variable = stress_xx
  []
  [syy]
    type = ElementAverageValue
    variable = stress_yy
  []
  [sxy]
    type = ElementAverageValue
    variable = stress_xy
  []
  [szz]
    type = ElementAverageValue
    variable = stress_zz
  []
  [sxz]
    type = ElementAverageValue
    variable = stress_xz
  []
  [syz]
    type = ElementAverageValue
    variable = stress_yz
  []
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-10
  start_time = 0
  end_time = ${t_end}
  [TimeStepper]
    type = TimeSequenceStepper
    time_sequence = '${t_levels}'
  []
[]
[Outputs]
  csv = true
  file_base = ${run_id}
[]
