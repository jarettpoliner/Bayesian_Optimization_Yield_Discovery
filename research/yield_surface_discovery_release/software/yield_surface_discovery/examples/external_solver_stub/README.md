# External-solver stub

`stub_solver.py` shows the contract between the adaptive controller and **any** command-line mechanics code:

* **Input:** the requested unit loading direction (Mandel z1, z2, z3; Cauchy components are also available as placeholders).
* **Output:** a JSON file with the first-yield stress (`stress_at_yield`, Cauchy s_xx, s_yy, s_xy) or `yield_radius`. Optional fields: `achieved_direction`, `status` and `metadata`.

Run the closed loop with the stub:

```
python scripts/run_closed_loop.py --config configs/example_solver.yaml
```

To use your own solver, change `oracle.command` in the config.

* Placeholders: `{z1} {z2} {z3} {s_xx} {s_yy} {s_xy} {run_dir} {run_id} {python}`.
* Each test runs in its own `run_dir`, with `stdout.log` and `stderr.log` kept.
* A non-zero exit code or a missing result marks that candidate unavailable, and the controller moves on to the next one.
