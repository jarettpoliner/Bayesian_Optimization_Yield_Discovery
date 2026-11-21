"""Writer r1 = proposed_inputs/moose_input_biphasic.py (unchanged copy) + one declared option, pin_mode:
  "origin"  : pin disp_x, disp_y at the node (0, 0) only  -- identical to psbench.moose_input (production configuration)
  "corners" : pin disp_x, disp_y at all four corner nodes (0,0), (L,0), (0,L), (L,L).
Rationale (STAGE_S1_VERIFICATION.md): the displacement variables are the PERIODIC FLUCTUATION (the macroscopic gradient enters through the
homogenization scalar), so the four corners are periodic images with identical fluctuation; pinning the origin already implies zero at all
corners.  On meshes where libMesh's periodic constraints leave one corner node untied (S2 h40: 160 instead of 162 constrained DOFs) the
"corners" pin restores the missing tie exactly; on correctly constrained meshes it is redundant."""
"""PROPOSED (not used by any campaign) MOOSE input writer for a plane-stress periodic RVE with TWO elastoplastic phases.

It is a block-for-block copy of psbench.moose_input.write_input (which is NOT modified), with these changes only:
  * each phase p in {matrix, inclusion} gets its own SolidMechanicsHardeningPowerRule (linear: value_0 = sy_p, epsilon0 = sy_p / H_p,
    exponent 1), its own SolidMechanicsPlasticJ2 and its own block-restricted ComputeMultiPlasticityStress wrapped by
    ComputeLagrangianWrappedStress (the existing matrix pattern, repeated on the inclusion block);
  * the eqps aux field (plastic_internal_parameter[0]) is defined on both blocks;
  * new postprocessors: eqps_i_avg/max, eqps_rve_avg (volume-weighted over the whole RVE), plastic volume fractions per phase
    (volume average of an indicator eqps > pthr), local plastic-work density per phase w_p = sy_p k + H_p k^2 / 2 (exact for associative J2
    with linear isotropic hardening and radial return, k = eqps), vm_i_max/avg, phase stress averages;
  * the Terminator stops on the chosen criterion variable (eqps_rve_avg by default).
Everything else (formulation, plane stress, periodic BCs, stress control, solver settings, outputs) is unchanged.
"""
import math
from pathlib import Path

SQRT2 = math.sqrt(2.0)
PHASES = ("matrix", "inclusion")


def _blk(name, typ, body=None, indent="  "):
    s = f"{indent}[{name}]\n{indent}  type = {typ}\n"
    for k, v in (body or {}).items():
        s += f"{indent}  {k} = {v}\n"
    return s + f"{indent}[]\n"


def write_input(path, mesh_file, spec, direction, times, stop_value, stop_pp="eqps_rve_avg", L=1.0, pthr=1e-6,
                field_output=False, checkpoint=False, restart_from=None, pin_mode="corners"):
    """spec: {'matrix': {E, nu, sy, H}, 'inclusion': {E, nu, sy, H}} -- both phases plastic."""
    dxx, dyy, dm = direction; dxy = dm / SQRT2
    PIN_COORDS = {"origin": "0 0 0", "corners": f"0 0 0; {L} 0 0; 0 {L} 0; {L} {L} 0"}[pin_mode]
    s = f"""# two-plastic-phase input (biphasic_benchmark, writer r1, pin_mode = {pin_mode}). Plane-stress periodic RVE, stress control.
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
    file = '{mesh_file}'
  []
  [pin]
    type = ExtraNodesetGenerator
    input = file
    coord = '{PIN_COORDS}'
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
    expression = '{float(dxx):.17g}*t'
  []
  [sxy_t]
    type = ParsedFunction
    expression = '{float(dxy):.17g}*t'
  []
  [syy_t]
    type = ParsedFunction
    expression = '{float(dyy):.17g}*t'
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
      translation = '{L} 0 0'
    []
    [y]
      variable = 'disp_x disp_y'
      primary = bottom
      secondary = top
      translation = '0 {L} 0'
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
"""
    tag = {"matrix": "m", "inclusion": "i"}
    s += _blk("strain", "ComputeLagrangianWPSStrain", dict(out_of_plane_strain="strain_zz", homogenization_gradient_names="homogenization_gradient"))
    s += _blk("hgrad", "ComputeHomogenizedLagrangianStrain")
    for ph in PHASES:
        p, t = spec[ph], tag[ph]
        s += _blk(f"C_{ph}", "ComputeIsotropicElasticityTensor", dict(block=ph, youngs_modulus=p["E"], poissons_ratio=p["nu"]))
        s += _blk(f"stress_{ph}", "ComputeLagrangianWrappedStress", dict(block=ph, rotate_old_stress="false"))
        s += _blk(f"stress_{ph}_base", "ComputeMultiPlasticityStress",
                  dict(block=ph, plastic_models=f"j2_{t}", ep_plastic_tolerance=1e-10, tangent_operator="nonlinear"))
    comps = dict(xx=(0, 0), yy=(1, 1), xy=(0, 1), zz=(2, 2), xz=(0, 2), yz=(1, 2))
    for c, (i, j) in comps.items():
        s += _blk(f"s{c}", "RankTwoCartesianComponent", dict(rank_two_tensor="cauchy_stress", property_name=f"s{c}", index_i=i, index_j=j))
    for c in ("xx", "yy", "xy", "zz"):
        i, j = comps[c]
        s += _blk(f"e{c}", "RankTwoCartesianComponent", dict(rank_two_tensor="total_strain", property_name=f"e{c}", index_i=i, index_j=j))
    s += _blk("vm", "RankTwoInvariant", dict(rank_two_tensor="cauchy_stress", property_name="vm", invariant="VonMisesStress"))
    s += "[]\n\n[Postprocessors]\n"
    ex = "'initial timestep_end'"; pp = ""
    for c in ("xx", "yy", "xy", "zz"):
        pp += _blk(f"s_{c}", "ElementAverageMaterialProperty", dict(mat_prop=f"s{c}", execute_on=ex))
        pp += _blk(f"e_{c}", "ElementAverageMaterialProperty", dict(mat_prop=f"e{c}", execute_on=ex))
    for c in ("zz", "xz", "yz"):
        for vt in ("max", "min"):
            pp += _blk(f"s{c}_qp_{vt}", "ElementExtremeMaterialProperty", dict(mat_prop=f"s{c}", value_type=vt, execute_on=ex))
    # volume-weighted averages of the element-projected eqps scalar field (ElementAverageValue = integral / volume)
    pp += _blk("eqps_rve_avg", "ElementAverageValue", dict(variable="eqps", execute_on=ex))
    pp += _blk("wp_rve_avg", "ElementAverageValue", dict(variable="wp", execute_on=ex))
    for ph in PHASES:
        t = tag[ph]
        pp += _blk(f"eqps_{t}_avg", "ElementAverageValue", dict(variable="eqps", block=ph, execute_on=ex))
        pp += _blk(f"eqps_{t}_max", "ElementExtremeValue", dict(variable="eqps", value_type="max", block=ph, execute_on=ex))
        pp += _blk(f"pvf_{t}", "ElementAverageValue", dict(variable="plastic_ind", block=ph, execute_on=ex))
        pp += _blk(f"wp_{t}_avg", "ElementAverageValue", dict(variable="wp", block=ph, execute_on=ex))
        pp += _blk(f"vm_{t}_max", "ElementExtremeMaterialProperty", dict(mat_prop="vm", value_type="max", block=ph, execute_on=ex))
        pp += _blk(f"vm_{t}_avg", "ElementAverageMaterialProperty", dict(mat_prop="vm", block=ph, execute_on=ex))
        for c in ("xx", "yy", "xy"):
            pp += _blk(f"s_{c}_{t}_avg", "ElementAverageMaterialProperty", dict(mat_prop=f"s{c}", block=ph, execute_on=ex))
        pp += _blk(f"vol_{t}", "VolumePostprocessor", dict(block=ph, execute_on="initial"))
    pp += _blk("vm_max", "ElementExtremeMaterialProperty", dict(mat_prop="vm", value_type="max", execute_on=ex))
    pp += _blk("nl_its", "NumNonlinearIterations", dict(execute_on="timestep_end"))
    pp += _blk("l_its", "NumLinearIterations", dict(execute_on="timestep_end"))
    s += pp + "[]\n\n[AuxVariables]\n"
    fields = ("vm", "sxx", "syy", "sxy", "szz") if field_output else ()
    for v in ("eqps", "plastic_ind", "wp") + tuple(f"{v}_f" for v in fields):
        s += f"  [{v}]\n    family = MONOMIAL\n    order = CONSTANT\n  []\n"
    s += "[]\n\n[AuxKernels]\n"
    s += _blk("eqps", "MaterialStdVectorAux", dict(variable="eqps", property="plastic_internal_parameter", index=0, block="'matrix inclusion'",
                                                    execute_on="'initial timestep_end'"))
    s += _blk("plastic_ind", "ParsedAux", dict(variable="plastic_ind", coupled_variables="eqps", expression=f"'if(eqps > {pthr:g}, 1, 0)'",
                                                execute_on="'initial timestep_end'"))
    for ph in PHASES:
        p, t = spec[ph], tag[ph]
        s += _blk(f"wp_{t}", "ParsedAux", dict(variable="wp", block=ph, coupled_variables="eqps",
                                                expression=f"'{float(p['sy']):.10g}*eqps + 0.5*{float(p['H']):.10g}*eqps*eqps'", execute_on="'initial timestep_end'"))
    for v in fields:
        s += _blk(f"{v}_f", "MaterialRealAux", dict(variable=f"{v}_f", property=v, execute_on="timestep_end"))
    s += "[]\n\n"
    if field_output:
        s += "[VectorPostprocessors]\n" + _blk("fields", "ElementValueSampler", dict(variable="'eqps vm_f sxx_f syy_f sxy_f szz_f'", sort_by="id", execute_on="final"))
        s += "[]\n\n"
    s += "[UserObjects]\n"
    for ph in PHASES:
        p, t = spec[ph], tag[ph]
        s += _blk(f"hardening_{t}", "SolidMechanicsHardeningPowerRule", dict(value_0=f"{float(p['sy']):.10g}", epsilon0=f"{float(p['sy'] / p['H']):.10g}", exponent=1))
        s += _blk(f"j2_{t}", "SolidMechanicsPlasticJ2", dict(yield_strength=f"hardening_{t}", yield_function_tolerance="1e-6", internal_constraint_tolerance="1e-10"))
    s += _blk("stop", "Terminator", dict(expression=f"'{stop_pp} > {float(stop_value):.6g}'", fail_mode="HARD", execute_on="timestep_end"))
    s += f"""[]

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
    time_sequence = '{' '.join(f'{float(t):.10g}' for t in times[1:])}'
  []
  start_time = 0
  end_time = {float(times[-1]):.10g}
[]

[Outputs]
  print_linear_residuals = false
  [csv]
    type = CSV
    file_base = '{Path(path).stem}_out'
    execute_vector_postprocessors_on = final
  []
"""
    if checkpoint:
        s += f"  [cp]\n    type = Checkpoint\n    num_files = {len(times) + 10}\n  []\n"
    s += "[]\n"
    if restart_from:
        s += f"\n[Problem]\n  restart_file_base = '{restart_from}'\n[]\n"
    Path(path).write_text(s)
