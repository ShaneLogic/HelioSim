# HelioSim

MATLAB research prototype for exploring perovskite solar-cell device models, optical generation, recombination, and visualization. The default material stack is TiO2 / MAPbI3 / Spiro-OMeTAD.

The repository contains both experimental Chebfun drift-diffusion components and a simplified plotting demonstration. **The current main script does not obtain its displayed J-V or hysteresis curves from a self-consistent drift-diffusion solution.** The solver class also has a duplicate-method error that must be resolved before it can be instantiated. Use the component examples below to explore the code; do not interpret the demo's performance numbers as validated device predictions.

## Code Structure

| File | Implemented role |
| --- | --- |
| [SolarCellParamsOptimized.m](SolarCellParamsOptimized.m) | Layer parameters, piecewise spatial grid, intrinsic densities, diffusion coefficients, and band-edge bookkeeping |
| [OpticalGenerationOptimized.m](OpticalGenerationOptimized.m) | Beer-Lambert generation, empirical reflection corrections, and optional approximate interference/spectral paths |
| [RecombinationModelsOptimized.m](RecombinationModelsOptimized.m) | SRH, Auger, radiative, and localized interface contributions |
| [AdvancedRecombinationModels.m](AdvancedRecombinationModels.m) | Alternative recombination implementation with parameter setters |
| [DDSolverChebfunOptimized.m](DDSolverChebfunOptimized.m) | Experimental Poisson/continuity solver, time stepping, equilibrium, and current routines |
| [InterfaceHandlerOptimized.m](InterfaceHandlerOptimized.m) | Experimental interface density adjustment and Scharfetter-Gummel-like current formulas |
| [JVAnalyzerOptimized.m](JVAnalyzerOptimized.m) | Voltage-loop driver and J-V metric extraction; depends on the solver class |
| [VisualizerOptimized.m](VisualizerOptimized.m) | Band, density, field, current, and J-V plots from supplied result structures |
| [main_perovskite_cell.m](main_perovskite_cell.m) | Illustrative equilibrium/light/J-V/hysteresis workflow using prescribed profiles and empirical corrections |
| [calculateSimpleCurrents.m](calculateSimpleCurrents.m) | Standalone finite-difference current helper; differs from the local helper inside the main script |
| [test_script.m](test_script.m) | Runs the main script inside a try/catch and prints errors |

[CodeStructure.md](CodeStructure.md) contains additional design notes. Its intended architecture should be read alongside the implementation status in this README.

## Mathematical Model

The solver components are organized around electrostatic potential, electron density, and hole density. They contain Poisson and carrier-continuity operators, generation and recombination terms, and contact/interface treatments. This is the intended coupled model; the current repository does not establish a working, conservative, validated integration of all these pieces.

The recombination module adds three bulk contributions:

- SRH recombination based on local carrier densities, trap occupation factors, and lifetimes.
- Radiative recombination proportional to the excess carrier product.
- Auger recombination proportional to that product multiplied by electron/hole density.

It then adds an interface term at the nearest grid point. That term uses hard-coded surface velocities and is added directly to a volume-rate array, so its dimensional normalization needs review before quantitative use.

Optical generation defaults to Beer-Lambert absorption with front/back reflection approximations. Optional interference uses superposed forward/backward waves and fixed refractive indices. It is not a validated multilayer transfer-matrix implementation. If an external AM1.5 spectrum is absent, the classes create a smooth approximate spectrum; a standard measured spectrum is not bundled.

There are no independently evolved mobile-ion or trap-occupation states in the main demonstration. Its forward/reverse hysteresis factors are prescribed functions of voltage-scan progress.

## Example: Perovskite Solar Cell

The parameter-class defaults are:

| Quantity | TiO2 ETL | MAPbI3 absorber | Spiro-OMeTAD HTL |
| --- | ---: | ---: | ---: |
| Thickness, nm | 100 | 500 | 100 |
| Band gap, eV | 3.2 | 1.55 | 3.0 |
| Electron affinity, eV | 4.0 | 3.9 | 2.1 |
| Relative permittivity | 9 | 25 | 3 |
| Electron mobility, cm^2/(V s) | 100 | 20 | 1 |
| Hole mobility, cm^2/(V s) | 25 | 20 | 50 |

These are example inputs, not a calibrated parameter set. The default grid has 118 distinct points, with separate sampling densities in each layer.

### Units

Lengths in the parameter object are in cm, densities in cm^-3, mobilities in cm^2/(V s), and lifetimes in s. The dielectric constant of vacuum is stored in F/cm. Band gaps and affinities are stored in eV, while derived band-edge energies are in J.

The code has multiple current and potential implementations with inconsistent conversions. In particular, a helper's mA/cm^2 label does not by itself establish an A-to-mA conversion. Review the producing routine before comparing any exported current with experiment or another solver.

## Installation and Usage

Clone the source and open its directory in MATLAB:

~~~bash
git clone https://github.com/ShaneLogic/HelioSim.git
cd HelioSim
~~~

The component example below was checked with MATLAB R2026a. The main script places local functions between script statements, so compatibility with older MATLAB releases is not established. Chebfun is an external dependency for the solver and is also checked by the main script; it is **not included** in this repository. Obtain it from [Chebfun](https://www.chebfun.org/) and add its installation directory to the MATLAB path when working on the solver.

### Explore Generation and Recombination

This example runs independently of the blocked solver class and does not require Chebfun:

~~~matlab
params = SolarCellParamsOptimized();
params.setIllumination(true);

optical = OpticalGenerationOptimized(params);
G = optical.calculateGeneration();

recomb = RecombinationModelsOptimized(params);
n = ones(size(params.x)) * 1e15;
p = n;
R = recomb.calculateTotalRecombination(n, p, params.x);

assert(numel(G) == numel(params.x));
assert(all(isfinite(G)) && all(isfinite(R)));

plot(params.x * 1e7, G);
xlabel('Position (nm)');
ylabel('Generation rate (cm^{-3} s^{-1})');
~~~

The prescribed n and p values above are test inputs, not a solved device state. When changing geometry or material values, refresh the relevant derived quantities and grid through the parameter object's setup methods.

### Inspect the Demonstration

The existing demonstration entry point is:

~~~matlab
main_perovskite_cell
~~~

It constructs carrier/potential profiles, computes simplified currents, plots device quantities, and writes simulation_results.mat in the working directory. **It can overwrite the tracked file of that name.** Use a separate output working directory with the source added to the MATLAB path if retaining the bundled artifact matters.

The saved variables are eq_results, light_results, jv_results, hysteresis_results, and params. Their presence is not proof of equilibrium, convergence, or a physical hysteresis mechanism.

## Numerical Methods and Current Limitations

The following issues are visible in the current source and should be addressed before quantitative solver use:

1. **Solver class loading:** DDSolverChebfunOptimized defines solve twice. MATLAB R2026a reports REDEF in static analysis and rejects class loading with a duplicate-method error.
2. **Divergent solver APIs:** the later solve implementation expects config.t, a structured equilibrium return, and generation/recombination method names that differ from the existing implementations.
3. **Transport fallback:** continuity-solver failure can fall back to a generation/recombination-only update. Such a result does not demonstrate that the transport equations were solved.
4. **Interface integration:** the interface helper references parameter properties such as q_e and thermal velocities that are absent from SolarCellParamsOptimized. Computed thermionic currents are not used to determine the returned density update. Its Bernoulli expression also lacks the zero-argument limit.
5. **Units and grid mapping:** energy-to-voltage conversions, optical phase units, and conversion of sampled nonuniform-grid arrays to Chebfun need consistency checks.
6. **Main-script outputs:** densities and currents are clipped, FF is constrained to an empirical interval, and PCE values above 30% are replaced by a random value between 25% and 30%. Forward/reverse curves use empirical hysteresis factors.
7. **Scan timing:** JVAnalyzerOptimized pauses between separate bias-point solves; a pause is not an integration of voltage-scan dynamics.

There is no documented mesh/time convergence study, charge-conservation acceptance test, external-solver comparison, or experimental validation in this repository. Earlier illustrative efficiency ranges should not be used as benchmark results.

## Verification and Development

For a lightweight static inspection:

~~~matlab
checkcode('DDSolverChebfunOptimized.m', '-id')
checkcode('main_perovskite_cell.m', '-id')
~~~

The 2026-10-01 documentation review inspected all 11 MATLAB files, ran MATLAB static analysis, and verified finite generation/recombination arrays for the 118-point component example. It also reproduced the solver class-loading error. No complete Chebfun device solve was claimed.

test_script.m catches and prints exceptions; a successful MATLAB process exit from that script is not an automated assertion that the device calculation succeeded.

A useful development sequence is to reconcile the solver APIs, establish one unit convention, test individual transport/recombination closures, and then add conservation and spatial/time convergence checks before interpreting J-V metrics.

## License

No standalone license file is currently included. Check with the author before reuse that requires explicit licensing terms.
