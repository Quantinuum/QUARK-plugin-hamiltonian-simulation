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

Benchmark for Hamiltonian dynamics for fermionic systems. Paper is `Phys. Rev. Research 7, 043146 (2025)`.

The system considered is a square lattice with width and height `lx` and `ly`, with open or periodic boundary conditions (determined by `periodic`), with one or two species of fermions (determined by `two_spin_species`). The Hamiltonian contains tight-binding hopping terms for all neighbours on the lattice. The initial state is a product state with a disk of sites at the center of the lattice filled with fermions, the other sites being empty. The system is evolved with the Hamiltonian for `n_trot` Trotter steps of size `dt`. One measure the fermion occupation number on all sites with `n_shots` shots, and computes the imbalance of fermions in/out the initial disk, namely the density of fermions within the initial disk minus outside the initial disk, and compare it to the exact value. The score output is the minimal number of two-qubit gates that an ideal perfect quantum computer would have to run to be able to detect an error in the output of the tested hardware by 3 standard deviations. `create_plot` determines whether plots of the imbalance are generated at the end of the benchmark.

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


Benchmark for Hamiltonian dynamics for fermionic systems and Fermionic Fourier Transform. Paper is `arXiv:2605.01440`.

The system considered is two coupled chains of size `l` with periodic boundary conditions, for one fermion species. The Hamiltonian contains tight-binding hopping terms for all sites of the first chain, and a coupling `epsilon` between the first chain and the second chain. The first chain is initialized in the ground state of the tigh-binding hopping terms with a Fermionic Fourier Transform (FFT). The first and second chains are then evolved for `n_trot` Trotter steps with size `dt`. A FFT is performed on the second chain and the density of fermions is measured in the second chain with `n_shots` shots. The output is related to the spectral function of the tight-binding terms on the first chain at a given value of frequency. The benchmark runs `l` different values of frequencies. Three scores are output, the fidelity between all measured and exact values, the fidelity with the dispersion relation obtained for noiseless, and the fidelity with the exact dispersion relation. Three modes are available for the compilation of the FFT, determined by `mode`, with `CZ_pyramid` involving only nearest-neighbour sites, `CX_ladder` minimizing the number of two-qubit gates, and `log2N` minimizing the two-qubit gate depth. `create_plot` determines whether plots of the spectral function are generated at the end of the benchmark.


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


Benchmark for state preparation and Hamiltonian dynamics for spin systems. Paper is `arXiv:2509.05206`.

The system considered is a square lattice with width and height `lx`, `ly`, and with periodic or open boundaries (determined by `periodic`). The Hamiltonian is an Ising model on this lattice, with XX coupling and Z magnetic field `h`. The benchmark performs an adiabatic evolution from the ground state of the Hamiltonian given by only X field on all sites, and the Ising model, using `n_trot` Trotter steps, and time steps that are scaled by a factor `dt`. The benchmark runs two different circuits, one where one starts from the ground state of the X field, and one where one starts from the first excited states with one qubit flipped, in order to compute a derivative of energy with respect to the initial entropy. Then the benchmark runs two additional mirror circuits that are used to estimate the amount of entropy injected by noise in the system. The benchmark outputs the energy of the state obtained, its entropy and its effective temperature calculated as explained in the reference. The circuits generated are all run with `n_shots` shots, so in total there are four times `n_shots` shots run. 


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
