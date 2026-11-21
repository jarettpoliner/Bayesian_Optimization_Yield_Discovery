# Analytical-oracle examples

* **Minimal example (about 10 s):** `python scripts/run_minimal_example.py --check`. This runs the physics-informed controller for 12 real tests against the Mg–4%Li analytical model, writes all outputs and a figure, and compares the result with `data/minimal_example/expected.json`.
* **Full analytical benchmarks:** `python scripts/run_benchmark.py --config configs/benchmark_physics_informed.yaml` (or `configs/benchmark_blind.yaml`).
* **Bit-exact reproduction of the registered runs:** use `configs/reproduce_registered_physics_informed.yaml` / `configs/reproduce_registered_blind.yaml`, which require the pinned environment.
* **Available analytical models** (`yield_surface_discovery.analytical.MODELS`):
  * `cb2004_mg4li`: the benchmark target;
  * `von_mises`: parameter `sigma_y`.
