# TL1 periodic truss RVE, stress-controlled proportional loading. direction (Mandel unit) = 1 0 0
[Mesh]
  displacements = 'disp_x disp_y'
  allow_renumbering = false
  [fmg]
    type = FileMeshGenerator
    file = '/home/jarettpoliner/projects/yield_surface_moose/plane_stress/closed_loop/truss_lattice_benchmark_2026/oracle/mesh_TL1/mesh.e'
    constraint_matrix = '/home/jarettpoliner/projects/yield_surface_moose/plane_stress/closed_loop/truss_lattice_benchmark_2026/oracle/mesh_TL1/constraint_matrix.m'
    allow_renumbering = false
  []
[]
[GlobalParams]
  displacements = 'disp_x disp_y'
[]
[Problem]
  # block 99 is the kernel-free dummy element that only hosts the control-node DOFs C1, C2
  material_coverage_check = false
[]
[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]
[AuxVariables]
  [area]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3 4 5 6'
  []
  [axial_stress]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3 4 5 6'
  []
  [plastic_stretch]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3 4 5 6'
  []
  [total_stretch]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3 4 5 6'
  []
  [hardening_variable]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3 4 5 6'
  []
[]
[AuxKernels]
  [area_1]
    type = ConstantAux
    variable = area
    block = 1
    value = 0.05
    execute_on = 'initial timestep_begin'
  []
  [area_2]
    type = ConstantAux
    variable = area
    block = 2
    value = 0.03
    execute_on = 'initial timestep_begin'
  []
  [area_3]
    type = ConstantAux
    variable = area
    block = 3
    value = 0.03
    execute_on = 'initial timestep_begin'
  []
  [area_4]
    type = ConstantAux
    variable = area
    block = 4
    value = 0.02
    execute_on = 'initial timestep_begin'
  []
  [area_5]
    type = ConstantAux
    variable = area
    block = 5
    value = 0.02
    execute_on = 'initial timestep_begin'
  []
  [area_6]
    type = ConstantAux
    variable = area
    block = 6
    value = 0.02
    execute_on = 'initial timestep_begin'
  []
  [axial_stress]
    type = MaterialRealAux
    property = axial_stress
    variable = axial_stress
    block = '1 2 3 4 5 6'
    execute_on = 'initial timestep_end'
  []
  [plastic_stretch]
    type = MaterialRealAux
    property = plastic_stretch
    variable = plastic_stretch
    block = '1 2 3 4 5 6'
    execute_on = 'initial timestep_end'
  []
  [total_stretch]
    type = MaterialRealAux
    property = total_stretch
    variable = total_stretch
    block = '1 2 3 4 5 6'
    execute_on = 'initial timestep_end'
  []
  [hardening_variable]
    type = MaterialRealAux
    property = hardening_variable
    variable = hardening_variable
    block = '1 2 3 4 5 6'
    execute_on = 'initial timestep_end'
  []
[]
[Kernels]
  [solid_disp_x]
    type = StressDivergenceTensorsTruss
    block = '1 2 3 4 5 6'
    component = 0
    variable = disp_x
    area = area
  []
  [solid_disp_y]
    type = StressDivergenceTensorsTruss
    block = '1 2 3 4 5 6'
    component = 1
    variable = disp_y
    area = area
  []
[]
[Materials]
  [mat_FRAME]
    type = PlasticTruss
    block = '1 2 3'
    youngs_modulus = 200000.0
    yield_stress = 300.0
    hardening_constant = 2000.0
    absolute_tolerance = 1e-12
    relative_tolerance = 1e-12
  []
  [mat_INFILL]
    type = PlasticTruss
    block = '4 5 6'
    youngs_modulus = 200000.0
    yield_stress = 200.0
    hardening_constant = 2000.0
    absolute_tolerance = 1e-12
    relative_tolerance = 1e-12
  []
[]
[BCs]
  [pin_x]
    type = DirichletBC
    variable = disp_x
    boundary = pin
    value = 0
  []
  [pin_y]
    type = DirichletBC
    variable = disp_y
    boundary = pin
    value = 0
  []
  [rot_gauge]
    type = DirichletBC
    variable = disp_y
    boundary = c1
    value = 0
  []
[]
[Functions]
  [f1x]
    type = ParsedFunction
    expression = '1.7320508075688772*t'
  []
  [f2x]
    type = ParsedFunction
    expression = '0.0*t'
  []
  [f2y]
    type = ParsedFunction
    expression = '0.0*t'
  []
  [dtfun]
    type = ParsedFunction
    expression = 'if(t < 7.678758572543162, 0.7678758580221922, 0.01616580753730931)'
  []
[]
[NodalKernels]
  [load_c1x]
    type = UserForcingFunctorNodalKernel
    variable = disp_x
    boundary = c1
    functor = f1x
  []
  [load_c2x]
    type = UserForcingFunctorNodalKernel
    variable = disp_x
    boundary = c2
    functor = f2x
  []
  [load_c2y]
    type = UserForcingFunctorNodalKernel
    variable = disp_y
    boundary = c2
    functor = f2y
  []
[]
[Postprocessors]
  [u0_x]
    type = NodalVariableValue
    variable = disp_x
    nodeid = 0
  []
  [u0_y]
    type = NodalVariableValue
    variable = disp_y
    nodeid = 0
  []
  [u1_x]
    type = NodalVariableValue
    variable = disp_x
    nodeid = 1
  []
  [u1_y]
    type = NodalVariableValue
    variable = disp_y
    nodeid = 1
  []
  [u2_x]
    type = NodalVariableValue
    variable = disp_x
    nodeid = 2
  []
  [u2_y]
    type = NodalVariableValue
    variable = disp_y
    nodeid = 2
  []
  [u3_x]
    type = NodalVariableValue
    variable = disp_x
    nodeid = 5
  []
  [u3_y]
    type = NodalVariableValue
    variable = disp_y
    nodeid = 5
  []
  [u4_x]
    type = NodalVariableValue
    variable = disp_x
    nodeid = 3
  []
  [u4_y]
    type = NodalVariableValue
    variable = disp_y
    nodeid = 3
  []
  [u5_x]
    type = NodalVariableValue
    variable = disp_x
    nodeid = 4
  []
  [u5_y]
    type = NodalVariableValue
    variable = disp_y
    nodeid = 4
  []
  [u6_x]
    type = NodalVariableValue
    variable = disp_x
    nodeid = 6
  []
  [u6_y]
    type = NodalVariableValue
    variable = disp_y
    nodeid = 6
  []
  [u7_x]
    type = NodalVariableValue
    variable = disp_x
    nodeid = 7
  []
  [u7_y]
    type = NodalVariableValue
    variable = disp_y
    nodeid = 7
  []
  [u8_x]
    type = NodalVariableValue
    variable = disp_x
    nodeid = 8
  []
  [u8_y]
    type = NodalVariableValue
    variable = disp_y
    nodeid = 8
  []
  [u9_x]
    type = NodalVariableValue
    variable = disp_x
    nodeid = 9
  []
  [u9_y]
    type = NodalVariableValue
    variable = disp_y
    nodeid = 9
  []
  [u10_x]
    type = NodalVariableValue
    variable = disp_x
    nodeid = 10
  []
  [u10_y]
    type = NodalVariableValue
    variable = disp_y
    nodeid = 10
  []
  [u11_x]
    type = NodalVariableValue
    variable = disp_x
    nodeid = 11
  []
  [u11_y]
    type = NodalVariableValue
    variable = disp_y
    nodeid = 11
  []
  [s0]
    type = ElementalVariableValue
    variable = axial_stress
    elementid = 0
  []
  [ep0]
    type = ElementalVariableValue
    variable = plastic_stretch
    elementid = 0
  []
  [eps0]
    type = ElementalVariableValue
    variable = total_stretch
    elementid = 0
  []
  [hv0]
    type = ElementalVariableValue
    variable = hardening_variable
    elementid = 0
  []
  [s1]
    type = ElementalVariableValue
    variable = axial_stress
    elementid = 1
  []
  [ep1]
    type = ElementalVariableValue
    variable = plastic_stretch
    elementid = 1
  []
  [eps1]
    type = ElementalVariableValue
    variable = total_stretch
    elementid = 1
  []
  [hv1]
    type = ElementalVariableValue
    variable = hardening_variable
    elementid = 1
  []
  [s2]
    type = ElementalVariableValue
    variable = axial_stress
    elementid = 2
  []
  [ep2]
    type = ElementalVariableValue
    variable = plastic_stretch
    elementid = 2
  []
  [eps2]
    type = ElementalVariableValue
    variable = total_stretch
    elementid = 2
  []
  [hv2]
    type = ElementalVariableValue
    variable = hardening_variable
    elementid = 2
  []
  [s3]
    type = ElementalVariableValue
    variable = axial_stress
    elementid = 3
  []
  [ep3]
    type = ElementalVariableValue
    variable = plastic_stretch
    elementid = 3
  []
  [eps3]
    type = ElementalVariableValue
    variable = total_stretch
    elementid = 3
  []
  [hv3]
    type = ElementalVariableValue
    variable = hardening_variable
    elementid = 3
  []
  [s4]
    type = ElementalVariableValue
    variable = axial_stress
    elementid = 4
  []
  [ep4]
    type = ElementalVariableValue
    variable = plastic_stretch
    elementid = 4
  []
  [eps4]
    type = ElementalVariableValue
    variable = total_stretch
    elementid = 4
  []
  [hv4]
    type = ElementalVariableValue
    variable = hardening_variable
    elementid = 4
  []
  [s5]
    type = ElementalVariableValue
    variable = axial_stress
    elementid = 5
  []
  [ep5]
    type = ElementalVariableValue
    variable = plastic_stretch
    elementid = 5
  []
  [eps5]
    type = ElementalVariableValue
    variable = total_stretch
    elementid = 5
  []
  [hv5]
    type = ElementalVariableValue
    variable = hardening_variable
    elementid = 5
  []
  [s6]
    type = ElementalVariableValue
    variable = axial_stress
    elementid = 6
  []
  [ep6]
    type = ElementalVariableValue
    variable = plastic_stretch
    elementid = 6
  []
  [eps6]
    type = ElementalVariableValue
    variable = total_stretch
    elementid = 6
  []
  [hv6]
    type = ElementalVariableValue
    variable = hardening_variable
    elementid = 6
  []
  [s7]
    type = ElementalVariableValue
    variable = axial_stress
    elementid = 7
  []
  [ep7]
    type = ElementalVariableValue
    variable = plastic_stretch
    elementid = 7
  []
  [eps7]
    type = ElementalVariableValue
    variable = total_stretch
    elementid = 7
  []
  [hv7]
    type = ElementalVariableValue
    variable = hardening_variable
    elementid = 7
  []
  [s8]
    type = ElementalVariableValue
    variable = axial_stress
    elementid = 8
  []
  [ep8]
    type = ElementalVariableValue
    variable = plastic_stretch
    elementid = 8
  []
  [eps8]
    type = ElementalVariableValue
    variable = total_stretch
    elementid = 8
  []
  [hv8]
    type = ElementalVariableValue
    variable = hardening_variable
    elementid = 8
  []
  [s9]
    type = ElementalVariableValue
    variable = axial_stress
    elementid = 9
  []
  [ep9]
    type = ElementalVariableValue
    variable = plastic_stretch
    elementid = 9
  []
  [eps9]
    type = ElementalVariableValue
    variable = total_stretch
    elementid = 9
  []
  [hv9]
    type = ElementalVariableValue
    variable = hardening_variable
    elementid = 9
  []
  [s10]
    type = ElementalVariableValue
    variable = axial_stress
    elementid = 10
  []
  [ep10]
    type = ElementalVariableValue
    variable = plastic_stretch
    elementid = 10
  []
  [eps10]
    type = ElementalVariableValue
    variable = total_stretch
    elementid = 10
  []
  [hv10]
    type = ElementalVariableValue
    variable = hardening_variable
    elementid = 10
  []
  [s11]
    type = ElementalVariableValue
    variable = axial_stress
    elementid = 11
  []
  [ep11]
    type = ElementalVariableValue
    variable = plastic_stretch
    elementid = 11
  []
  [eps11]
    type = ElementalVariableValue
    variable = total_stretch
    elementid = 11
  []
  [hv11]
    type = ElementalVariableValue
    variable = hardening_variable
    elementid = 11
  []
  [E11]
    type = ParsedPostprocessor
    pp_names = 'u10_x u11_x u11_y'
    expression = '(u10_x/2.0)'
    use_t = true
  []
  [E22]
    type = ParsedPostprocessor
    pp_names = 'u10_x u11_x u11_y'
    expression = '(u11_y/1.7320508075688772)'
    use_t = true
  []
  [E12m]
    type = ParsedPostprocessor
    pp_names = 'u10_x u11_x u11_y'
    expression = '(1.4142135623730951*0.5*((u11_x-1.0*u10_x/2.0)/1.7320508075688772))'
    use_t = true
  []
  [ep_norm]
    type = ParsedPostprocessor
    pp_names = 'u10_x u11_x u11_y'
    expression = 'sqrt(((u10_x/2.0)-t*(0.0001237179148263*1.0+-4.1230690636828e-05*0.0+-2.8571432257468e-08*0.0))^2+((u11_y/1.7320508075688772)-t*(-4.1242706443397e-05*1.0+0.0001676901405822*0.0+-3.8095278801901e-09*0.0))^2+((1.4142135623730951*0.5*((u11_x-1.0*u10_x/2.0)/1.7320508075688772))-t*(-6.098637220231e-20*1.0+1.6940658945086e-20*0.0+0.0002309401301109*0.0))^2)'
    use_t = true
  []
  [epL_norm]
    type = ParsedPostprocessor
    pp_names = 'ep0 ep1 ep2 ep3 ep4 ep5 ep6 ep7 ep8 ep9 ep10 ep11'
    expression = 'sqrt((0.35714285714286187*ep0+-3.8459253727671054e-13*ep1+-3.8459253727671054e-13*ep2+0.35714285714286187*ep3+-2.5639502485114036e-13*ep4+-2.5639502485114036e-13*ep5+0.14285714285714474*ep6+-3.8459253727671054e-13*ep7+-2.5639502485114036e-13*ep8+0.14285714285714474*ep9+-2.5639502485114036e-13*ep10+-3.8459253727671054e-13*ep11)^2+(-0.11902275169007538*ep0+0.19998955652890124*ep1+0.19998955652890124*ep2+-0.11902275169007538*ep3+0.1333263710192675*ep4+0.1333263710192675*ep5+-0.04760910067603015*ep6+0.19998955652890124*ep7+0.1333263710192675*ep8+-0.04760910067603015*ep9+0.1333263710192675*ep10+0.19998955652890124*ep11)^2+(-8.247862041711447e-05*ep0+0.24498362719794592*ep1+-0.2449143451558744*ep2+-8.247862041711447e-05*ep3+0.16332241813196394*ep4+-0.16327623010391626*ep5+-3.299144816684579e-05*ep6+0.24498362719794592*ep7+-0.16327623010391626*ep8+-3.299144816684579e-05*ep9+0.16332241813196394*ep10+-0.2449143451558744*ep11)^2)'
  []
[]
[UserObjects]
  [stop]
    type = Terminator
    expression = 'epL_norm > 0.00015000000000000001'
    fail_mode = HARD
    error_level = NONE
    execute_on = 'timestep_end'
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
  petsc_options_value = 'lu petsc'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-09
  nl_max_its = 30
  l_tol = 1e-12
  l_max_its = 50
  end_time = 404.1451884327327
  [TimeStepper]
    type = FunctionDT
    function = dtfun
    cutback_factor_at_failure = 0.5
  []
[]
[Outputs]
  csv = true
  exodus = false
  [pc]
    type = PerfGraphOutput
    level = 0
    execute_on = final
  []
[]
