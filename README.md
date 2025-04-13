# HelioSim: Drift-Diffusion Solar Cell Simulator

A comprehensive MATLAB toolkit for simulating the physics and performance of multi-layer solar cells using the drift-diffusion model.

![Solar Cell Simulation](https://i.imgur.com/example.png)

## Features

- Detailed drift-diffusion solver for solar cell operation
- Multi-layer device support (ETL/Absorber/HTL architecture)
- Comprehensive physical models for carrier generation and recombination
- Band diagram visualization
- J-V curve analysis and solar cell performance metrics
- Parameter variation for device optimization

## Physical Model

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

### Physical Processes

The simulator includes detailed models for:

- **Optical Generation**: Beer-Lambert absorption model with optional AM1.5G spectrum
- **Recombination Mechanisms**:
  - Shockley-Read-Hall (SRH) recombination
  - Radiative (band-to-band) recombination
  - Auger recombination
- **Interface Physics**: Band offsets and interface recombination
- **Band Diagram Calculation**: Including band bending at interfaces

## Installation

1. Clone this repository or download all files
2. Open MATLAB (version 2016b or newer recommended)
3. Navigate to the HelioSim directory
4. Run the main script to start the simulation:
   ```matlab
   runSimulation
   ```

## Quick Start Guide

### Running Your First Simulation

1. **Basic Equilibrium Simulation**:
   ```matlab
   % Create parameters with default values
   params = SolarCellParams();
   
   % Configure simulation settings
   config = struct('t_start', 0, 't_end', 1e-9, 'num_time_steps', 51, ...
                   'rel_tol', 1e-6, 'abs_tol', 1e-8, ...
                   'illumination', false, 'voltage_sweep', false);
   
   % Create solver and run simulation
   solver = DDSolver(params, config);
   results = solver.solve();
   
   % Visualize results
   visualizer = Visualizer(params);
   visualizer.plotBandDiagram(results);
   ```

2. **Generate J-V Curve**:
   ```matlab
   % Set illumination to true
   params.setIllumination(true);
   
   % Create J-V analyzer and generate curve
   analyzer = JVAnalyzer(params, solver);
   jv_results = analyzer.generateJVCurve(-0.2, 1.2, 20);
   
   % Display performance metrics
   fprintf('Voc: %.4f V\n', jv_results.Voc);
   fprintf('Jsc: %.4f mA/cm^2\n', jv_results.Jsc);
   fprintf('Fill Factor: %.4f\n', jv_results.FF);
   fprintf('Efficiency: %.2f%%\n', jv_results.PCE);
   ```

### Using the GUI

For easier interaction, use the built-in GUI by running:
```matlab
runSimulation
```

This opens an interactive interface where you can:
- Set device parameters
- Choose simulation type (equilibrium, applied bias, illumination, J-V)
- Visualize results
- Perform parameter sweeps for optimization

## Customizing Device Parameters

You can customize your solar cell device by modifying parameters:

```matlab
% Create default parameters
params = SolarCellParams();

% Modify absorber layer properties
params.setAbsorberThickness(500e-7);  % 500 nm
params.setAbsorberBandgap(1.55);      % 1.55 eV
params.setAbsorberMobility(10);       % 10 cm²/Vs

% Modify contact layers
params.setETLDoping(5e17);            % Donor concentration in ETL
params.setHTLDoping(5e17);            % Acceptor concentration in HTL

% Update illumination state
params.setIllumination(true);

% Set applied voltage
params.setAppliedVoltage(0.5);        % 0.5 V
```

## Code Structure

- `DDSolver.m`: Core drift-diffusion solver
- `SolarCellParams.m`: Device parameter definitions
- `OpticalGeneration.m`: Light absorption and carrier generation
- `RecombinationModels.m`: Carrier recombination mechanisms
- `InterfaceHandler.m`: Interface condition handling
- `Visualizer.m`: Result visualization
- `JVAnalyzer.m`: J-V curve generation and performance analysis
- `runSimulation.m`: Main script with optional GUI
- `main_solar_cell.m`: Simple command-line example

## Example: Parameter Sweep

```matlab
% Create parameter and solver objects
params = SolarCellParams();
config = struct('t_start', 0, 't_end', 1e-9, 'num_time_steps', 51, ...
                'rel_tol', 1e-6, 'abs_tol', 1e-8, ...
                'illumination', true, 'voltage_sweep', false);
solver = DDSolver(params, config);
analyzer = JVAnalyzer(params, solver);

% Define bandgap range
bandgaps = 1.1:0.1:1.7;  % 1.1 to 1.7 eV
PCE = zeros(size(bandgaps));

% Sweep over bandgaps
for i = 1:length(bandgaps)
    params.setAbsorberBandgap(bandgaps(i));
    jv_results = analyzer.generateJVCurve();
    PCE(i) = jv_results.PCE;
    fprintf('Bandgap = %.2f eV, PCE = %.2f%%\n', bandgaps(i), PCE(i));
end

% Plot results
figure;
plot(bandgaps, PCE, 'o-', 'LineWidth', 2);
xlabel('Bandgap (eV)');
ylabel('Power Conversion Efficiency (%)');
title('PCE vs. Absorber Bandgap');
grid on;
```

## Advanced Usage

### Implementing Custom Generation Profiles

You can modify the optical generation profile by extending the `OpticalGeneration` class:

```matlab
% Create a custom generation profile
optical_gen = OpticalGeneration(params);
custom_gen = @(x) params.G_max * exp(-(x-params.x(1)).^2/1e-10);
G = arrayfun(custom_gen, params.x);

% Use in simulation
% [your simulation code]
```

### Using Different Recombination Models

You can modify recombination parameters to model different materials:

```matlab
% Create recombination model
recomb = RecombinationModels(params);

% Modify parameters for perovskite-like material
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