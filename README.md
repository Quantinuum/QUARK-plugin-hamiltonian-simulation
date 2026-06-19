# QUARK-plugin-hamiltonian-simulation

This is a QUARK plugin that implements hamiltonian simulation benchmarks. It should be used within the [QUARK framework](https://github.com/QUARK-framework/QUARK-framework)

See the directory `config_examples` for usage examples.

Run a benchmark config using:

```terminal
quark -c <config_file>
```

## Benchmarks

The descriptions of the benchmarks implemented within this QUARK plugin and their configuration APIs are given below.

Example benchmark configurations can be found in [config_examples](config_examples/).

### Free Fermion Benchmark

Pipeline step name: `free_fermion_benchmark`

Defined in [free_fermion_benchmark.py](src/quark_plugin_hamiltonian_simulation/free_fermion_benchmark/free_fermion_benchmark.py).

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `lx` | `int` | `2` | Number of lattice sites in the x direction. |
| `ly` | `int` | `2` | Number of lattice sites in the y direction. |
| `periodic` | `bool` | `False` | Whether to use periodic boundary conditions. |
| `two_spin_species` | `bool` | `False` | Whether to include two spin species in the model. |
| `create_plot` | `bool` | `False` | Whether to display result plots and stabilizer histograms during postprocessing. |
| `n_trot` | `int` | `2` | Number of Trotter steps. |
| `n_shots` | `int` | `200` | Number of shots used for each generated circuit. |
| `dt` | `float` | `0.5` | Trotter time step. |

Example:

```yaml
- "free_fermion_benchmark": {lx: 4, ly: 4, periodic: false, n_trot: 4, two_spin_species: false, create_plot: true}
```

### Spectral Function Benchmark

Pipeline step name: `spectral_function_benchmark`

Defined in [spectral_function.py](src/quark_plugin_hamiltonian_simulation/spectral_function_benchmark/spectral_function.py).

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `l` | `int` | `8` | System size parameter used by the benchmark circuits. |
| `dt` | `float` | `0.5` | Trotter time step. |
| `n_trot` | `int` | `10` | Number of Trotter steps. |
| `n_shots` | `int` | `200` | Number of shots used for each generated circuit. |
| `mode` | `str` | `"log2N"` | Circuit construction mode. Accepted values are `"log2N"`, `"CZ_pyramid"`, and `"CX_ladder"`. |
| `create_plot` | `bool` | `False` | Whether to display the measured and exact spectral-function plots during postprocessing. |
| `epsilon` | `float` | `np.pi * 0.8 / (n_trot * dt)` | Coupling parameter used to build the circuit sequence. |

Example:

```yaml
- "spectral_function_benchmark": {l: 8, n_trot: 10, dt: 0.5, epsilon: 0.5, n_shots: 100, mode: "log2N", create_plot: true}
```

### Thermal Adiabatic Square Ising Benchmark

Pipeline step name: `thermal_adiabatic_square_ising`

Defined in [thermal_adiabatic_square_ising.py](src/quark_plugin_hamiltonian_simulation/thermal_adiabatic_benchmark/thermal_adiabatic_square_ising.py).

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `lx` | `int` | `2` | Number of lattice sites in the x direction. |
| `ly` | `int` | `2` | Number of lattice sites in the y direction. |
| `h` | `float` | `1` | Transverse-field parameter. |
| `periodic` | `bool` | `False` | Whether to use periodic boundary conditions. |
| `dt` | `float` | `0.2` | Trotter time step. |
| `n_trot` | `int` | `2` | Number of Trotter steps. |
| `n_shots` | `int` | `200` | Number of shots used for each generated circuit. |

Example:

```yaml
- "thermal_adiabatic_square_ising": {lx: 4, ly: 4, periodic: false, n_trot: 10, h: 2, dt: 0.2}
```
