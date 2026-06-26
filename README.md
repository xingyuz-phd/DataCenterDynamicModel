# Dynamic Modeling of Data-Center Power Delivery for Power-System Resonance Analysis

This repository contains the supplementary MATLAB code for the manuscript:

> Xingyu Zhao and Junbo Zhao, **"Dynamic Modeling of Data-Center Power Delivery for Power System Resonance Analysis,"** arXiv:2604.06624v1, 2026.  
> Manuscript: <https://arxiv.org/html/2604.06624v1>

The code implements the numerical studies associated with the paper's component-informed data-center power-delivery model and its integration into positive-sequence power-system dynamic analysis. The repository is code-only: generated figures, generated CSV outputs, MATLAB autosave files, and raw workload traces were intentionally removed.

## Paper context

The paper develops a dynamic model of the online, double-conversion data-center power-delivery chain and uses it for power-system oscillation studies. The modeled chain includes:

1. an active front-end (AFE) rectifier connected to the grid-side point of common coupling (PCC),
2. a UPS DC-link capacitor,
3. a voltage-source inverter (VSI),
4. an aggregated PSU-array equivalent, and
5. a downstream DC-DC/load equivalent representing server-side CPU/GPU demand.

The model is formulated as a time-invariant positive-sequence representation so that it can be combined with phasor-domain grid models and reduced to a small-signal state-space model. The scripts in this repository support the main numerical workflows in the paper:

- equilibrium computation of the nonlinear data-center/grid model,
- bandwidth-based PI-controller tuning,
- comparison between full three-phase and reduced positive-sequence/QSS models,
- eigenvalue and participation-factor analysis,
- power oscillation amplification (POA) analysis from server-load disturbance to grid-visible active power, and
- time-domain simulations driven by sinusoidal or GPU-workload load variations.

## Repository layout

```text
TPWRS_DataCenter_SupplementalCode/
├── README.md
├── .gitignore
├── .gitattributes
└── src/
    ├── case9/
    │   ├── case9.m
    │   ├── datacenter_init_from_pf_qss.m
    │   ├── datacenter_port_model.m
    │   ├── default_datacenter_params.m
    │   ├── gfl_init_from_pf_qss_newton.m
    │   ├── gfl_port_model.m
    │   ├── gfm_init_from_pf_qss_newton.m
    │   ├── gfm_port_model.m
    │   ├── sm_init_from_pf_fsolve.m
    │   └── sm_port_model.m
    └── sdcib/
        ├── SDCIB.m
        ├── SDCIB_simplified_eigs_step.m
        ├── ModelValidation.m
        └── gpu_load_response.m
```

## Code-to-paper map

| Paper component | Main script(s) | Description |
| --- | --- | --- |
| Component-level data-center model | `src/sdcib/SDCIB.m`, `src/case9/datacenter_port_model.m` | Implements the reduced positive-sequence data-center power-delivery model with AFE, DC link, VSI, PSU, and downstream DC-DC/load states. |
| Model validation in the SDCIB case | `src/sdcib/ModelValidation.m` | Compares the full cascaded three-phase representation with the reduced QSS/positive-sequence model under a load step. This corresponds to the validation study in the SDCIB section of the paper. |
| SDCIB modal, participation, and POA analysis | `src/sdcib/SDCIB.m` | Computes the equilibrium, linearizes the model, evaluates eigenvalues and participation factors, computes POA curves, and performs parameter scans over controller bandwidth, load level, and grid strength. |
| Simplified SDCIB load-step example | `src/sdcib/SDCIB_simplified_eigs_step.m` | Provides a compact eigenvalue/participation-factor workflow and a reduced load-step response example. |
| Realistic GPU-load time-domain response in SDCIB | `src/sdcib/gpu_load_response.m` | Uses an external GPU-load trace to simulate workload-induced propagation from server load to PCC power and to compute FFT spectra. |
| Modified 3-machine 9-bus grid-integration case | `src/case9/case9.m` | Builds the modified 9-bus system with a synchronous machine, a grid-forming inverter, a grid-following inverter, and the data-center load at Bus 8; then performs modal/POA analysis and GPU-load time-domain simulation. |
| SM/GFM/GFL/data-center dynamic port models | `src/case9/*_port_model.m`, `src/case9/*_init_*.m` | Initialization and dynamic port-equation routines used by the 9-bus study. |

## Main workflows

### 1. Single-data-center infinite-bus case

The SDCIB scripts study the intrinsic oscillatory behavior of the data-center power-delivery chain when the data center is connected to an infinite bus through a Thevenin impedance. This workflow is useful for isolating data-center internal modes before coupling the model to a larger grid.

Recommended run order:

```matlab
cd src/sdcib
ModelValidation
SDCIB_simplified_eigs_step
SDCIB
```

`ModelValidation.m` compares the full three-phase model and the reduced QSS positive-sequence model. `SDCIB.m` then computes eigenvalues, participation factors, POA curves, controller/load/grid-strength sensitivity scans, and sinusoidal load-response examples.

### 2. Realistic GPU-workload response

Two scripts use a time-varying GPU-load trace:

```matlab
cd src/sdcib
gpu_load_response
```

```matlab
cd src/case9
case9
```

These scripts require an external file named `GPU_data.csv` in the current MATLAB working directory. The expected columns are:

```text
t_seconds,P_gpu_W
```

The raw workload trace is not included in this code-only release. When using your own trace, keep the same column names or modify the corresponding `readtable` section in the scripts.

### 3. Modified 3-machine 9-bus case

The 9-bus workflow is implemented in:

```matlab
cd src/case9
case9
```

This script:

1. solves the power flow of the modified 9-bus system,
2. performs Kron reduction for the retained dynamic ports,
3. initializes the synchronous-machine, GFM, GFL, and data-center port models,
4. forms the reduced small-signal model,
5. computes eigenvalues and modal participation/coupling information,
6. computes POA transfer gains from data-center load variation to multiple grid-side ports, and
7. runs a GPU-load-driven time-domain simulation when `GPU_data.csv` is available.

In the included setup, the data-center load is connected at Bus 8. The three dynamic generation units are represented by a synchronous-machine model, a grid-forming inverter model, and a grid-following inverter model.

## MATLAB requirements

Recommended environment:

- MATLAB R2022b or newer.
- Optimization Toolbox, required for `fsolve` and `optimoptions`.

The scripts also use standard MATLAB functionality including `ode15s`, `readtable`, `writetable`, `eig`, `fft`, plotting routines, and script-local functions. No Simulink model or third-party package is included.

This repository was prepared as a clean code release. The numerical cases were not re-run during packaging.

## Generated files

The following files may be generated when running the scripts and are ignored by Git:

```text
eigvals_sorted.csv
participation_abs_norm.csv
*.fig
*.png
*.pdf
*.asv
```

If exact figure reproduction is needed, run the scripts from their respective folders and save the generated figures manually.

## Notes on modifying experiments

Common parameters are defined near the beginning of each main script. Useful entries include:

- `p.p_load0`, `p.p_load1`, and `p.sin_amp` in the SDCIB scripts,
- controller bandwidth targets in the `tg` structure of `SDCIB.m`,
- grid impedance or SCR-related parameters in `SDCIB.m`,
- `scale_lev` and Bus 8 data-center loading in `case9.m`, and
- the parameter structure returned by `src/case9/default_datacenter_params.m`.

For a new workload trace, place the data in `GPU_data.csv` or adjust the file name and column names in the scripts that call `readtable`.

## Citation

If this code is used in academic work, please cite the associated manuscript:

```bibtex
@article{zhao2026datacenter,
  title   = {Dynamic Modeling of Data-Center Power Delivery for Power System Resonance Analysis},
  author  = {Zhao, Xingyu and Zhao, Junbo},
  journal = {arXiv preprint arXiv:2604.06624},
  year    = {2026}
}
```

## License

No license file is included in this packaged code release. Add a `LICENSE` file before making the repository public if a specific open-source license is intended.
