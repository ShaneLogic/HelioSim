%% main_solar_cell.m
% Main script for running the multi-layer solar cell simulation
% Handles ETL, Absorber, and HTL layers with proper interfaces

clear; close all; clc;

disp('Initializing multi-layer solar cell simulation...');

% Load or create solar cell parameters
params = SolarCellParams();

% Configure simulation settings
sim_config = struct(...
    't_start', 0, ...
    't_end', 1e-9, ...
    'num_time_steps', 51, ...
    'rel_tol', 1e-6, ...
    'abs_tol', 1e-8, ...
    'illumination', true, ...
    'voltage_sweep', true);

% Create solver instance
solver = DDSolver(params, sim_config);

% Run equilibrium simulation
disp('Solving drift-diffusion equations for equilibrium state...');
eq_results = solver.solve();

% Visualize equilibrium results
visualizer = Visualizer(params);
visualizer.plotSteadyState(eq_results);
visualizer.plotBandDiagram(eq_results);

% Run illuminated simulation if configured
if sim_config.illumination
    disp('Running simulation with illumination...');
    params.setIllumination(true);
    light_results = solver.solve();
    visualizer.plotComparisonWithLight(eq_results, light_results);
end

% Run J-V curve analysis if configured
if sim_config.voltage_sweep
    disp('Generating J-V characteristics...');
    analyzer = JVAnalyzer(params, solver);
    jv_results = analyzer.generateJVCurve();
    
    % Display performance metrics
    fprintf('\nSolar Cell Performance Metrics:\n');
    fprintf('Open-circuit voltage (Voc): %.4f V\n', jv_results.Voc);
    fprintf('Short-circuit current (Jsc): %.4f mA/cm^2\n', jv_results.Jsc);
    fprintf('Fill Factor (FF): %.4f\n', jv_results.FF);
    fprintf('Power Conversion Efficiency (PCE): %.2f%%\n', jv_results.PCE);
end

disp('Simulation completed.');