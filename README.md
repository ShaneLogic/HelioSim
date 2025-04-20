# HelioSim: Drift-Diffusion Solar Cell Simulator for Perovskite Devices

A comprehensive MATLAB toolkit for simulating the physics and performance of perovskite solar cells using the drift-diffusion model with Chebfun spectral methods.

## Table of Contents

1. [Code Structure](#code-structure)
2. [Mathematical Model](#mathematical-model)
3. [Numerical Methods](#numerical-methods)
4. [Example: Perovskite Solar Cell](#example-perovskite-solar-cell)
5. [Installation and Usage](#installation-and-usage)
6. [Advanced Features](#advanced-features)

## Code Structure

HelioSim is organized into several core modules that work together to simulate solar cell behavior:

### Core Components

1. **SolarCellParamsOptimized.m**
   - Contains physical constants, material parameters, and grid settings
   - Specialized for perovskite solar cells
   - Provides convenient parameter setting methods

2. **DDSolverChebfunOptimized.m**
   - Drift-diffusion equation solver using Chebfun spectral methods
   - Implements high-precision solutions for Poisson and continuity equations
   - Integrates interface handling and Scharfetter-Gummel discretization
   - Uses adaptive time stepping and Newton iteration

3. **RecombinationModelsOptimized.m**
   - Implements SRH, Auger, and radiative recombination models
   - Special handling for interface recombination
   - Supports energy-dependent trap models

4. **OpticalGenerationOptimized.m**
   - Carrier generation model using Beer-Lambert absorption
   - Supports wavelength-dependent absorption coefficients
   - Optimized for perovskite materials

5. **JVAnalyzerOptimized.m**
   - J-V curve analyzer
   - Calculates Voc, Jsc, fill factor, and efficiency
   - Implements maximum power point tracking
   - Optimized voltage scanning algorithm

6. **VisualizerOptimized.m**
   - Results visualization tool
   - Band diagram plotting
   - Carrier density visualization
   - J-V curve plotting
   - Electric field distribution visualization

7. **main_perovskite_cell.m**
   - Main script for perovskite solar cell simulation
   - Sets parameters and configuration
   - Runs simulation
   - Analyzes and visualizes results

### Data Flow

1. `main_perovskite_cell.m` creates a `SolarCellParamsOptimized` instance and sets parameters
2. Creates a `DDSolverChebfunOptimized` instance with parameters and configuration
3. Solver internally creates `RecombinationModelsOptimized` and `OpticalGenerationOptimized` instances
4. Solver solves drift-diffusion equations and returns results
5. Creates `JVAnalyzerOptimized` instance to analyze performance
6. Creates `VisualizerOptimized` instance to visualize results

## Mathematical Model

HelioSim implements a comprehensive drift-diffusion model that solves the coupled semiconductor equations:

### Core Equations

1. **Poisson's Equation**: 
   ```
   ∇²φ = q/ε₀ε * (p - n + N⁺ - N⁻)
   ```
   where φ is the electrostatic potential, q is the elementary charge, ε is the relative permittivity, n and p are electron and hole concentrations, and N⁺/N⁻ are ionized donor/acceptor concentrations.

2. **Carrier Continuity Equations**:
   ```
   ∂n/∂t = ∇·(Dn∇n + μn·n·∇φ) + G - R
   ∂p/∂t = ∇·(Dp∇p - μp·p·∇φ) + G - R
   ```
   where D is the diffusion coefficient, μ is the carrier mobility, G is the generation rate, and R is the recombination rate.

3. **Current Densities**:
   ```
   Jn = q·μn·n·E + q·Dn·∇n
   Jp = q·μp·p·E - q·Dp·∇p
   ```
   where E is the electric field.

### Boundary Conditions

- **Ohmic Contacts**: Fixed carrier concentrations at the electrodes
- **Interface Conditions**: Continuity of electric displacement and quasi-Fermi levels at material interfaces

### Physical Processes

The simulator includes detailed models for:

- **Optical Generation**: Beer-Lambert absorption model with AM1.5G spectrum
- **Recombination Mechanisms**:
  - Shockley-Read-Hall (SRH) recombination
  - Radiative (band-to-band) recombination
  - Auger recombination
  - Interface recombination
- **Interface Physics**: Band offsets and interface recombination
- **Band Diagram Calculation**: Including band bending at interfaces

## Numerical Methods

HelioSim employs advanced numerical techniques for accurate and efficient simulation:

### Spatial Discretization

- **Chebfun Spectral Methods**: Replaces finite differences with spectral methods
- **Adaptive Grid Refinement**: At interfaces for better resolution
- **High-Precision Derivative Calculation**: For accurate field and current calculations

### Time Integration

- **Implicit Time Stepping**: For numerical stability
- **Adaptive Time Step Control**: Based on solution dynamics
- **Newton Method**: For solving nonlinear equation systems

### Interface Handling

- **Scharfetter-Gummel Discretization**: For accurate current calculation
- **Band Discontinuity Treatment**: Accounts for energy band offsets
- **Thermionic Emission Model**: For carrier transport across interfaces

### Recombination Models

- **Position-Dependent Parameters**: For different material regions
- **Interface-Specific Recombination**: Enhanced recombination at interfaces
- **Trap Energy Distribution**: For realistic defect modeling

## Example: Perovskite Solar Cell

### Device Structure

The default simulation models a typical perovskite solar cell with the following structure:

- **ETL**: TiO2 (100 nm)
- **Absorber**: MAPbI3 (500 nm)
- **HTL**: Spiro-OMeTAD (100 nm)

### Key Parameters

#### TiO2 (ETL)
- Bandgap: 3.2 eV
- Electron affinity: 4.0 eV
- Dielectric constant: 9.0
- Electron mobility: 100 cm²/Vs
- Donor doping: 1×10¹⁷ cm⁻³

#### MAPbI3 (Absorber)
- Bandgap: 1.55 eV
- Electron affinity: 3.9 eV
- Dielectric constant: 25.0
- Electron/hole mobility: 20 cm²/Vs
- Intrinsic (undoped)
- Carrier lifetime: 1 μs

#### Spiro-OMeTAD (HTL)
- Bandgap: 3.0 eV
- Electron affinity: 2.1 eV
- Dielectric constant: 3.0
- Hole mobility: 50 cm²/Vs
- Acceptor doping: 1×10¹⁷ cm⁻³

#### Interface Parameters
- ETL/Absorber interface recombination velocity: 10⁴ cm/s
- Absorber/HTL interface recombination velocity: 10⁴ cm/s

### Simulation Results

The simulation produces the following key performance metrics for the perovskite solar cell:

- **Open-circuit voltage (Voc)**: ~1.0-1.1 V
- **Short-circuit current density (Jsc)**: ~22-24 mA/cm²
- **Fill factor (FF)**: ~0.75-0.80
- **Power conversion efficiency (PCE)**: ~18-22%

The simulation also generates detailed visualizations including:

- Band diagram showing band bending at interfaces
- Carrier concentration profiles
- Electric field distribution
- J-V characteristic curve
- Recombination rate profiles

## Installation and Usage

### Requirements

- MATLAB (version 2016b or newer recommended)
- Chebfun library (included in the package)

### Installation

1. Clone this repository or download all files
2. Open MATLAB
3. Navigate to the HelioSim directory
4. Run the main script to start the simulation:
   ```matlab
   main_perovskite_cell
   ```

### Basic Usage

```matlab
% Run the main perovskite cell simulation
main_perovskite_cell

% To modify parameters before running:
params = SolarCellParamsOptimized();
params.L_absorber = 500e-7;  % 500 nm absorber thickness
params.Eg_abs = 1.55;        % 1.55 eV bandgap
params.mu_n_abs = 20;        % 20 cm²/Vs electron mobility
params.mu_p_abs = 20;        % 20 cm²/Vs hole mobility

% Configure simulation
config.t_max = 1e-9;          % Maximum simulation time
config.illumination = true;   % Enable illumination

% Create solver and run
solver = DDSolverChebfunOptimized(params, config);
results = solver.solve();

% Analyze J-V characteristics
analyzer = JVAnalyzerOptimized(params, solver);
jv_results = analyzer.generateJVCurve(-0.1, 1.2, 30);

% Visualize results
visualizer = VisualizerOptimized(params);
visualizer.setResults(results);
visualizer.plotBandDiagram();
visualizer.setJVResults(jv_results);
visualizer.plotJVCurve();
```

## Advanced Features

### Parameter Sweeps

```matlab
% Example: Sweep absorber thickness
thicknesses = [300, 400, 500, 600, 700] * 1e-7;  % nm to cm
PCE = zeros(size(thicknesses));

for i = 1:length(thicknesses)
    params.L_absorber = thicknesses(i);
    % Update grid
    params.generateGrid();
    % Run simulation and get J-V results
    jv_results = analyzer.generateJVCurve();
    PCE(i) = jv_results.PCE;
    fprintf('Thickness = %.0f nm, PCE = %.2f%%\n', thicknesses(i)*1e7, PCE(i));
end

% Plot results
figure;
plot(thicknesses*1e7, PCE, 'o-', 'LineWidth', 2);
xlabel('Absorber Thickness (nm)');
ylabel('Power Conversion Efficiency (%)');
title('PCE vs. Absorber Thickness');
grid on;
```

### Hysteresis Analysis

The simulator can model hysteresis effects in perovskite solar cells by performing forward and reverse voltage scans:

```matlab
% Configure hysteresis analysis
config.scan_rate = 100;  % V/s
config.preconditioning = true;

% Run forward and reverse scans
[forward_results, reverse_results] = analyzer.performHysteresisAnalysis();

% Calculate hysteresis index
HI = analyzer.calculateHysteresisIndex(forward_results, reverse_results);
fprintf('Hysteresis Index: %.4f\n', HI);

% Visualize hysteresis
visualizer.plotHysteresisJVCurve(forward_results, reverse_results);
```

### Custom Optical Generation

```matlab
% Create custom generation profile
optical_gen = OpticalGenerationOptimized(params);
optical_gen.enable_interference = true;  % Enable interference effects
custom_gen = optical_gen.calculateGeneration();

% Use in simulation
solver.setGenerationProfile(custom_gen);
```

### Interface Engineering

```matlab
% Modify interface properties
params.S_ETL_abs = 1e3;  % ETL/absorber interface recombination velocity (cm/s)
params.S_abs_HTL = 1e3;  % absorber/HTL interface recombination velocity (cm/s)

% Add interface dipole
params.dipole_ETL_abs = 0.1;  % 0.1 eV dipole at ETL/absorber interface
recomb.B_rad_abs = 5e-10;        % Stronger radiative recombination
recomb.Cn_auger_abs = 1e-29;     % Faster Auger recombination
recomb.Cp_auger_abs = 1e-29;
```

## Contributing

Contributions to HelioSim are welcome! Please feel free to submit pull requests or create issues for bugs and feature requests.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- The drift-diffusion model implementation is based on semiconductor physics principles described in [Semiconductor Device Physics and Design by Umesh K. Mishra and Jasprit Singh](https://link.springer.com/book/10.1007/978-1-4020-6481-4)
- The solar cell architecture is inspired by typical perovskite and thin-film solar cell designs

## Contact

For questions and support, please open an issue on the GitHub repository page. 