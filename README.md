# Supplementary MATLAB Code for arXiv:2604.06624v1

This repository contains the MATLAB source code extracted from the supplementary archive for the paper at <https://arxiv.org/html/2604.06624v1>. The package is code-only: generated figures, generated CSV files, backup files, and raw data traces were intentionally removed.

## Repository layout

```text
TPWRS_DataCenter_SupplementalCode/
├── README.md
├── .gitignore
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

## Code-path relationships

| Path | Role | Main dependencies |
| --- | --- | --- |
| `src/case9/case9.m` | IEEE 9-bus workflow: power flow, Kron reduction, SM/GFM/GFL/data-center initialization, modal analysis, POA gain, and GPU-load time-domain simulation. | Calls the device initialization/model files in `src/case9/`. The GPU-load section requires `GPU_data.csv` in the working directory. |
| `src/case9/datacenter_*.m` | Data-center port model and quasi-steady-state initialization. | Used by `case9.m`. |
| `src/case9/gfl_*.m` | Grid-following inverter initialization and dynamic port model. | Used by `case9.m`. |
| `src/case9/gfm_*.m` | Grid-forming inverter initialization and dynamic port model. | Used by `case9.m`. |
| `src/case9/sm_*.m` | Synchronous-machine initialization and dynamic port model. | Used by `case9.m`. |
| `src/sdcib/SDCIB.m` | Single-data-center infinite-bus modal, POA, parameter-scan, and sinusoidal-response study. | Self-contained script with local functions. Writes eigenvalue/participation CSV outputs in the working directory. |
| `src/sdcib/SDCIB_simplified_eigs_step.m` | Reduced SDCIB eigenvalue, participation-factor, and load-step example. | Self-contained script with local functions. Writes eigenvalue/participation CSV outputs in the working directory. |
| `src/sdcib/ModelValidation.m` | Validation script comparing the full abc cascaded model and reduced QSS uv representation. | Self-contained script with local functions. |
| `src/sdcib/gpu_load_response.m` | Time-domain simulation and FFT analysis driven by a measured GPU load trace. | Self-contained script with local functions. Requires `GPU_data.csv` in the working directory. |

## MATLAB requirements

- Recommended MATLAB release: R2022b or newer.
- Required toolbox: Optimization Toolbox, for `fsolve` and `optimoptions`.
- Standard MATLAB functions used include `ode15s`, `readtable`, `writetable`, `eig`, plotting utilities, and local functions in script files.
- No Simulink models or third-party packages are included in this code-only release.

The code was organized and packaged without re-running the numerical cases in this environment.

## Data and generated outputs

The uploaded archive contained PDFs, CSV files, and a MATLAB autosave file. These items were removed to keep this repository code-only.

Two scripts require a GPU-load trace when running their time-domain sections:

- `src/case9/case9.m`
- `src/sdcib/gpu_load_response.m`

Place a file named `GPU_data.csv` in the same directory from which the script is executed. The expected columns are:

```text
t_seconds,P_gpu_W
```

`SDCIB.m` and `SDCIB_simplified_eigs_step.m` may generate `eigvals_sorted.csv` and `participation_abs_norm.csv`; these are treated as generated outputs and are excluded by `.gitignore`.

## Running the scripts

From MATLAB, either start in the target script folder or add the source tree to the path:

```matlab
addpath(genpath('src'));
```

Example runs:

```matlab
cd src/case9
case9
```

```matlab
cd src/sdcib
SDCIB
SDCIB_simplified_eigs_step
ModelValidation
% gpu_load_response requires GPU_data.csv in this folder.
gpu_load_response
```

## Notes for GitHub upload

- `.gitignore` excludes generated figures, CSV outputs, MATLAB autosave files, and temporary files.
- Add a license file before public distribution if the code will be released under a specific license.
- When citing this code, cite the associated paper at <https://arxiv.org/html/2604.06624v1>.
