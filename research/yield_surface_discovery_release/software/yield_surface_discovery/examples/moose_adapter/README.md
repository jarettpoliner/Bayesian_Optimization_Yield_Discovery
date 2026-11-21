# MOOSE adapter (optional)

> **Status: NOT EXECUTED WITH MOOSE.** The adapter (`src/yield_surface_discovery/solvers/moose.py`) is unit-tested against a fake executable that writes MOOSE-format CSV output (`tests/test_moose_adapter.py`). The input template `single_element_yield.i` has never been run by MOOSE. Check its object names against your MOOSE revision before use.

**One virtual experiment** is a homogeneous single element loaded proportionally under stress control along the requested direction:

* tractions on the four lateral faces; free ±z faces (exact plane stress);
* the load factor λ equals the Mandel radius of the applied stress, so the yield radius is read directly from it.

**First yield** is the first λ at which the element-averaged equivalent plastic strain (`eqps` column) exceeds `eqps_threshold` (1e-8):

1. a 20-step ramp to `lambda_cap`, with one ×4 capacity escalation;
2. bisection to a relative bracket ≤ 1e-4.

The achieved stress direction and the out-of-plane residuals are checked on every test. Every MOOSE run keeps its input, CSV, stdout and stderr under `work_dir/test_####/`.

**Hidden material.** `material_parameters` in the config are substituted into the template on the oracle side only. The adaptive controller receives only:

* requested direction;
* achieved direction;
* stress at yield;
* yield radius;
* status;
* run metadata.

**To run on a machine with MOOSE:**

1. Build or obtain a MOOSE application containing the Solid Mechanics module (for example the official container `idaholab/moose`, or a source build).
2. Set `oracle.executable` in `configs/example_moose.yaml`. If MOOSE runs in a container, set `launcher` to your container wrapper.
3. Check the template with `<app> --check-input -i examples/moose_adapter/single_element_yield.i` after filling the placeholders by hand once.
4. Run `python scripts/run_closed_loop.py --config configs/example_moose.yaml`.

**First validation target:** J2 plasticity. The analytical plane-stress von Mises radius is the cross-check (`yield_surface_discovery.analytical.von_mises_radius`).
