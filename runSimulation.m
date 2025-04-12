%% runSimulation.m
% Convenience script to run the solar cell simulation with various options

%% Choose simulation type
simType = menu('Select simulation type:', ...
               'Equilibrium', ...
               'Applied bias', ...
               'Illumination', ...
               'J-V Curve', ...
               'Parameter variation', ...
               'Exit');

switch simType
    case 1 % Equilibrium
        % Load parameters
        params = SolarCellParams();
        
        % Configure simulation
        config = struct('t_start', 0, ...
                        't_end', 1e-9, ...
                        'num_time_steps', 51, ...
                        'rel_tol', 1e-6, ...
                        'abs_tol', 1e-8, ...
                        'illumination', false, ...
                        'voltage_sweep', false);
        
        % Set voltage to equilibrium (0V)
        params.setAppliedVoltage(0);
        
        % Create solver and run simulation
        solver = DDSolver(params, config);
        results = solver.solve();
        
        % Visualize results
        visualizer = Visualizer(params);
        visualizer.plotSteadyState(results);
        visualizer.plotBandDiagram(results);
        
    case 2 % Applied bias
        % Load parameters
        params = SolarCellParams();
        
        % Get voltage from user
        voltage = inputdlg('Enter applied voltage (V):', 'Applied Bias', 1, {'0.5'});
        voltage = str2double(voltage{1});
        
        % Configure simulation
        config = struct('t_start', 0, ...
                        't_end', 1e-9, ...
                        'num_time_steps', 51, ...
                        'rel_tol', 1e-6, ...
                        'abs_tol', 1e-8, ...
                        'illumination', false, ...
                        'voltage_sweep', false);
        
        % Set applied voltage
        params.setAppliedVoltage(voltage);
        
        % Create solver and run simulation
        solver = DDSolver(params, config);
        results = solver.solve();
        
        % Visualize results
        visualizer = Visualizer(params);
        visualizer.plotSteadyState(results);
        visualizer.plotBandDiagram(results);
        visualizer.plotElectricField(results);
        visualizer.plotCurrentDensities(results);
        
    case 3 % Illumination
        % Load parameters
        params = SolarCellParams();
        
        % Configure simulation
        config = struct('t_start', 0, ...
                        't_end', 1e-9, ...
                        'num_time_steps', 51, ...
                        'rel_tol', 1e-6, ...
                        'abs_tol', 1e-8, ...
                        'illumination', true, ...
                        'voltage_sweep', false);
        
        % Set applied voltage
        voltage = inputdlg('Enter applied voltage (V):', 'Applied Bias', 1, {'0'});
        voltage = str2double(voltage{1});
        params.setAppliedVoltage(voltage);
        
        % Create solver
        solver = DDSolver(params, config);
        
        % Run dark simulation first
        params.setIllumination(false);
        dark_results = solver.solve();
        
        % Run illuminated simulation
        params.setIllumination(true);
        light_results = solver.solve();
        
        % Visualize results
        visualizer = Visualizer(params);
        visualizer.plotComparisonWithLight(dark_results, light_results);
        
    case 4 % J-V Curve
        % Load parameters
        params = SolarCellParams();
        
        % Configure simulation
        config = struct('t_start', 0, ...
                        't_end', 1e-9, ...
                        'num_time_steps', 51, ...
                        'rel_tol', 1e-6, ...
                        'abs_tol', 1e-8, ...
                        'illumination', true, ...
                        'voltage_sweep', true);
        
        % Create solver
        solver = DDSolver(params, config);
        
        % Create J-V analyzer
        jv_analyzer = JVAnalyzer(params, solver);
        
        % Get voltage range from user
        prompt = {'Start voltage (V):', 'End voltage (V):', 'Number of points:'};
        defaults = {'-0.2', '1.2', '20'};
        answer = inputdlg(prompt, 'J-V Curve Settings', 1, defaults);
        
        if ~isempty(answer)
            V_start = str2double(answer{1});
            V_end = str2double(answer{2});
            V_points = round(str2double(answer{3}));
            
            % Generate J-V curve
            jv_results = jv_analyzer.generateJVCurve(V_start, V_end, V_points);
        end
        
    case 5 % Parameter variation
        % Load parameters
        params = SolarCellParams();
        
        % Choose parameter to vary
        paramType = menu('Select parameter to vary:', ...
                         'Absorber thickness', ...
                         'Absorber bandgap', ...
                         'ETL doping', ...
                         'HTL doping', ...
                         'Mobility', ...
                         'Back to main menu');
        
        if paramType == 6
            runSimulation;
            return;
        end
        
        % Configure simulation
        config = struct('t_start', 0, ...
                        't_end', 1e-9, ...
                        'num_time_steps', 51, ...
                        'rel_tol', 1e-6, ...
                        'abs_tol', 1e-8, ...
                        'illumination', true, ...
                        'voltage_sweep', false);
        
        % Get number of variations
        n_var = inputdlg('Number of variations to simulate:', 'Variation Count', 1, {'3'});
        n_var = round(str2double(n_var{1}));
        
        % Initialize arrays for results
        Voc_array = zeros(1, n_var);
        Jsc_array = zeros(1, n_var);
        FF_array = zeros(1, n_var);
        PCE_array = zeros(1, n_var);
        param_values = zeros(1, n_var);
        
        % Create solver
        solver = DDSolver(params, config);
        
        % Create J-V analyzer
        jv_analyzer = JVAnalyzer(params, solver);
        
        % Loop through variations
        for i = 1:n_var
            % Set parameter value based on type
            switch paramType
                case 1 % Absorber thickness
                    % Values from 200 nm to 700 nm
                    thickness = 2e-5 + (i-1) * 5e-5/(n_var-1);
                    params.L_absorber = thickness;
                    param_values(i) = thickness * 1e4; % in µm
                    param_name = 'Absorber Thickness (µm)';
                    
                case 2 % Absorber bandgap
                    % Values from 1.2 eV to 2.0 eV
                    bandgap = 1.2 + (i-1) * 0.8/(n_var-1);
                    params.Eg_abs = bandgap;
                    params.ni_abs = sqrt(params.Nc_abs * params.Nv_abs) * ...
                                    exp(-bandgap/(2*(params.kb*params.T/params.q)));
                    param_values(i) = bandgap;
                    param_name = 'Absorber Bandgap (eV)';
                    
                case 3 % ETL doping
                    % Values from 1e16 to 1e18 cm^-3
                    doping = 1e16 * 10^((i-1) * 2/(n_var-1));
                    params.Nd_ETL = doping;
                    param_values(i) = log10(doping);
                    param_name = 'ETL Doping log_{10}(cm^{-3})';
                    
                case 4 % HTL doping
                    % Values from 1e16 to 1e18 cm^-3
                    doping = 1e16 * 10^((i-1) * 2/(n_var-1));
                    params.Na_HTL = doping;
                    param_values(i) = log10(doping);
                    param_name = 'HTL Doping log_{10}(cm^{-3})';
                    
                case 5 % Mobility
                    % Values from 0.1 to 10 cm^2/Vs
                    mobility = 0.1 * 10^((i-1) * 2/(n_var-1));
                    params.mu_n_abs = mobility;
                    params.mu_p_abs = mobility;
                    params.D_n(params.idx_absorber) = mobility * params.Vt;
                    params.D_p(params.idx_absorber) = mobility * params.Vt;
                    param_values(i) = mobility;
                    param_name = 'Absorber Mobility (cm^2/Vs)';
            end
            
            % Reinitialize grid if geometry changed
            if paramType == 1
                params.setupGrid();
                params.setupMaterialParameters();
            end
            
            % Set up solver
            solver = DDSolver(params, config);
            jv_analyzer = JVAnalyzer(params, solver);
            
            % Generate J-V curve
            jv_results = jv_analyzer.generateJVCurve(-0.2, 1.2, 10);
            
            % Store results
            Voc_array(i) = jv_results.Voc;
            Jsc_array(i) = abs(jv_results.Jsc);
            FF_array(i) = jv_results.FF;
            PCE_array(i) = jv_results.PCE;
        end
        
        % Plot parameter variation results
        figure('Name', 'Parameter Variation Results', 'Position', [100, 100, 1200, 800]);
        
        % Plot Voc
        subplot(2,2,1);
        plot(param_values, Voc_array, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
        xlabel(param_name);
        ylabel('V_{oc} (V)');
        title('Effect on Open-Circuit Voltage');
        grid on;
        
        % Plot Jsc
        subplot(2,2,2);
        plot(param_values, Jsc_array, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
        xlabel(param_name);
        ylabel('J_{sc} (mA/cm^2)');
        title('Effect on Short-Circuit Current');
        grid on;
        
        % Plot FF
        subplot(2,2,3);
        plot(param_values, FF_array, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
        xlabel(param_name);
        ylabel('Fill Factor');
        title('Effect on Fill Factor');
        grid on;
        
        % Plot PCE
        subplot(2,2,4);
        plot(param_values, PCE_array, 'ko-', 'LineWidth', 2, 'MarkerSize', 8);
        xlabel(param_name);
        ylabel('PCE (%)');
        title('Effect on Power Conversion Efficiency');
        grid on;
        
    case 6 % Exit
        return;
end