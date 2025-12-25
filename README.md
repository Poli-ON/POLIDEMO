# POLIDEMO: hh

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.apenergy.2025.126744-blue)](https://doi.org/10.1016/j.apenergy.2025.126744)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2023%2B-orange)](https://www.mathworks.com/products/matlab.html)

POLIDEMO is an open-source physics-based degradation model implemented in MATLAB for predicting the evolution of electrical and mechanical performance of lithium-ion batteries during usage. The framework combines electrochemical and mechanical modeling approaches to accurately predict battery degradation, including capacity fade and internal resistance increase.

## Table of Contents

- [Key Features](#key-features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Extending POLIDEMO](#extending-polidemo)
- [Citing POLIDEMO](#citing-polidemo)
- [Contributing](#contributing)
- [Contact](#contact)

## Key Features

The main contributions of POLIDEMO are:

- **Accurate knee point prediction** in the capacity loss curve, by modeling the battery internal resistance as a function of SEI (Solid Electrolyte Interphase) thickness and Loss of Active Material (LAM) in both electrodes.
- **Novel mechanical degradation law** quantifying LAM in the electrodes.
- **Efficient fatigue modeling** - Alternative to Paris' law to compute the increase in total crack area caused by mechanical fatigue, reducing both computational cost and the number of model parameters.
- **Physics-based irreversible swelling** modeling, improving the accuracy of parameter estimation by providing an additional observable for model identification.
- **Model parameters estimation** is improved by fitting not just the capacity decay curves, but also the degradation indicators and the irreversible swelling.
- **Cycle-jump acceleration technique** integration to accelerate lifetime simulations, particularly beneficial during parameter estimation routines that require multiple consecutive model evaluations.

## Requirements

- POLIDEMO requires **MATLAB** R2023 or newer. Earlier MATLAB releases have not been tested and may not be fully compatible
- No additional toolboxes required (uses only built-in MATLAB functions)

## Installation

You can obtain POLIDEMO in two ways:

### Option 1 — Clone the repository

```bash
git clone https://github.com/Poli-ON/POLIDEMO.git
cd POLIDEMO
```

### Option 2 — Download as ZIP (no Git required)

1. Download the latest ZIP package from the [releases page](https://github.com/Poli-ON/POLIDEMO/releases) or click [here](https://github.com/Poli-ON/POLIDEMO/archive/refs/heads/main.zip)
2. Extract the ZIP file to your desired location

### Add to MATLAB Path

After obtaining POLIDEMO, open MATLAB and add it to your path:

```matlab
% Navigate to the POLIDEMO directory
cd('path/to/POLIDEMO')

% Add POLIDEMO and all subdirectories to MATLAB path
addpath(genpath(pwd))
```

To make this permanent, you can save the path:

```matlab
savepath
```

## Quick Start

1. Clone or download POLIDEMO
2. Add POLIDEMO to your MATLAB path (see [Installation](#installation))
3. Navigate to the `Examples` folder
4. Run the main script:

```matlab
POLIDEMO_main
```

## Usage

### Reproducing Publication Results

To reproduce the results presented in the associated publication ([DOI: 10.1016/j.apenergy.2025.126744](https://doi.org/10.1016/j.apenergy.2025.126744)):

1. Open MATLAB
2. Navigate to the `Examples` folder
3. Execute the main script:
   ```matlab
   POLIDEMO_main
   ```

### Simulation Workflow

When you run `POLIDEMO_main.m`, the script automatically performs the following steps:

1. **Parameter Initialization**
   - Loads all model parameters via `Parameters_init.m`
   - Sets up electrochemical and mechanical model configurations

2. **Simulation Execution**
   - Runs the electrochemical–mechanical degradation simulation via `run_n_cycles.m`
   - For each cycle, `spm_simulation.m` is called, which:
     - Initializes electrochemical and mechanical state variables (`spm_initial_conditions.m`)
     - Integrates the Single Particle Model (SPM) ODEs using MATLAB's `ode15s` solver
     - Processes and saves simulation results (`spm_post_processing.m`)

3. **Acceleration**
   - Applies the cycle-jump acceleration technique via `spm_cycle_jump.m` to speed up lifetime simulations

### Output

The simulation generates results including:
- Capacity fade over cycles
- Internal resistance evolution
- SEI thickness growth
- Loss of Active Material (LAM) in both electrodes
- Irreversible swelling measurements

## Extending POLIDEMO

The model is designed to be modular and easy to extend. Here are some ways to customize POLIDEMO:

### Adding New Degradation Mechanisms

- Add new degradation mechanisms by creating additional functions inside `spm_odefun_aging_new.m`
- Follow the existing function structure and integrate with the main ODE solver

### Modifying Parameters

- Modify parameters and initial conditions in `Parameters_init.m`
- Adjust electrochemical, mechanical, or aging parameters as needed for your specific battery chemistry

### Code Structure

- Explore commented sections within each function to understand the code structure and customize behavior
- The modular design allows for easy modification of individual components without affecting the entire framework

### ⚠️ Important Notes

**Experimental Features:** Some advanced features are planned for future releases and must **not** be activated in the current version, as doing so will cause the code to fail:

- `EnableTemperature` – to consider thermal effects on aging
- `EnableCEI` – for CEI (Cathode Electrolyte Interphase) formation modeling
- `EnablePlating` – for lithium plating modeling

These features will be fully integrated and stable in upcoming updates of POLIDEMO.

## Citing POLIDEMO

If you use POLIDEMO in your research, please cite the following article:

**BibTeX:**
```bibtex
@article{Pistorio2025,
  title={POLIDEMO: An electrochemical-mechanical framework for modeling lithium-ion batteries degradation},
  author={Pistorio, Francesca and Clerici, Davide and Somà, Aurelio},
  journal={Applied Energy},
  year={2025},
  volume={XXX},
  pages={XXX--XXX},
  doi={10.1016/j.apenergy.2025.126744},
  url={https://github.com/Poli-ON/POLIDEMO}
}
```

**Plain Text:**
Pistorio, F., Clerici, D., & Somà, A. (2025). POLIDEMO: An electrochemical-mechanical framework for modeling lithium-ion batteries degradation. *Applied Energy*, XXX, XXX-XXX. https://doi.org/10.1016/j.apenergy.2025.126744

## Contributing

We welcome contributions to POLIDEMO! Here's how you can help:

- **Bug Reports**: Found a bug? Please open an issue describing the problem, steps to reproduce, and your MATLAB version.
- **Feature Requests**: Have an idea for a new feature? Open an issue to discuss it.
- **Code Contributions**: Fork the repository, make your changes, and submit a pull request. Please ensure your code follows the existing style and includes appropriate comments.

## License

Please refer to the LICENSE file in this repository for license information.

## Contact

For questions, comments, suggestions, or bug reports, please contact the maintainers:

- **Francesca Pistorio** - [francesca.pistorio@polito.it](mailto:francesca.pistorio@polito.it)
- **Davide Clerici** - [davide.clerici@polito.it](mailto:davide.clerici@polito.it)

## Acknowledgments

POLIDEMO is developed at Politecnico di Torino. We thank all contributors and users for their feedback and support.

