classdef Visualizer < handle
    % Visualizer - Class to create plots from simulation results
    % Creates various visualizations for solar cell analysis
    
    properties
        params      % Solar cell parameters
        lineStyles  % Line styles for different regions
        colors      % Colors for different regions
    end
    
    methods
        function obj = Visualizer(params)
            % Constructor - initialize with solar cell parameters
            obj.params = params;
            
            % Define line styles and colors for different regions
            obj.lineStyles = {'-', '-', '-'};
            obj.colors = {'b', 'r', 'g'};  % ETL, Absorber, HTL
        end
        
        function plotSteadyState(obj, results)
            % Plot steady-state potential, carrier densities, and currents
            
            % Extract data from results
            x = obj.params.x;
            phi = results.phi(end,:)';
            n = results.n(end,:)';
            p = results.p(end,:)';
            Jn = results.Jn;
            Jp = results.Jp;
            J_total = results.J_total;
            
            % Get indices for different regions
            idx_ETL = obj.params.idx_ETL;
            idx_abs = obj.params.idx_absorber;
            idx_HTL = obj.params.idx_HTL;
            
            % Create figure for potential and carrier densities
            figure('Name', 'Steady State Results', 'Position', [100, 100, 1200, 800]);
            
            % Plot electric potential
            subplot(2, 2, 1);
            hold on;
            plot(x(idx_ETL)*1e4, phi(idx_ETL), 'LineWidth', 2, 'Color', obj.colors{1});
            plot(x(idx_abs)*1e4, phi(idx_abs), 'LineWidth', 2, 'Color', obj.colors{2});
            plot(x(idx_HTL)*1e4, phi(idx_HTL), 'LineWidth', 2, 'Color', obj.colors{3});
            hold off;
            title('Electric Potential');
            xlabel('Position (µm)');
            ylabel('Potential (V)');
            legend('ETL', 'Absorber', 'HTL', 'Location', 'best');
            grid on;
            
            % Plot carrier densities (log scale)
            subplot(2, 2, 2);
            hold on;
            semilogy(x(idx_ETL)*1e4, n(idx_ETL), 'LineWidth', 2, 'Color', obj.colors{1}, 'LineStyle', '-');
            semilogy(x(idx_abs)*1e4, n(idx_abs), 'LineWidth', 2, 'Color', obj.colors{2}, 'LineStyle', '-');
            semilogy(x(idx_HTL)*1e4, n(idx_HTL), 'LineWidth', 2, 'Color', obj.colors{3}, 'LineStyle', '-');
            
            semilogy(x(idx_ETL)*1e4, p(idx_ETL), 'LineWidth', 2, 'Color', obj.colors{1}, 'LineStyle', '--');
            semilogy(x(idx_abs)*1e4, p(idx_abs), 'LineWidth', 2, 'Color', obj.colors{2}, 'LineStyle', '--');
            semilogy(x(idx_HTL)*1e4, p(idx_HTL), 'LineWidth', 2, 'Color', obj.colors{3}, 'LineStyle', '--');
            hold off;
            title('Carrier Densities');
            xlabel('Position (µm)');
            ylabel('Carrier Density (cm^{-3})');
            legend('ETL n', 'Absorber n', 'HTL n', 'ETL p', 'Absorber p', 'HTL p', 'Location', 'best');
            grid on;
            ylim([1e2, 1e18]);
            
            % Plot current densities
            subplot(2, 2, 3);
            hold on;
            plot(x*1e4, Jn*1e3, 'b-', 'LineWidth', 2);
            plot(x*1e4, Jp*1e3, 'r-', 'LineWidth', 2);
            plot(x*1e4, J_total*1e3, 'k-', 'LineWidth', 2);
            hold off;
            title('Current Densities');
            xlabel('Position (µm)');
            ylabel('Current Density (mA/cm^2)');
            legend('Electron', 'Hole', 'Total', 'Location', 'best');
            grid on;
            
            % Plot electric field
            subplot(2, 2, 4);
            hold on;
            plot(x(idx_ETL)*1e4, results.E_field(idx_ETL)/1e3, 'LineWidth', 2, 'Color', obj.colors{1});
            plot(x(idx_abs)*1e4, results.E_field(idx_abs)/1e3, 'LineWidth', 2, 'Color', obj.colors{2});
            plot(x(idx_HTL)*1e4, results.E_field(idx_HTL)/1e3, 'LineWidth', 2, 'Color', obj.colors{3});
            hold off;
            title('Electric Field');
            xlabel('Position (µm)');
            ylabel('Electric Field (kV/cm)');
            legend('ETL', 'Absorber', 'HTL', 'Location', 'best');
            grid on;
        end
        
        function plotBandDiagram(obj, results)
            % Plot energy band diagram
            
            % Extract data
            x = obj.params.x;
            Ec = results.Ec;
            Ev = results.Ev;
            Efn = results.Efn;
            Efp = results.Efp;
            
            % Get indices for different regions
            idx_ETL = obj.params.idx_ETL;
            idx_abs = obj.params.idx_absorber;
            idx_HTL = obj.params.idx_HTL;
            
            % Convert to eV for plotting
            Ec_eV = Ec / obj.params.q;
            Ev_eV = Ev / obj.params.q;
            Efn_eV = Efn / obj.params.q;
            Efp_eV = Efp / obj.params.q;
            
            % Create figure
            figure('Name', 'Band Diagram', 'Position', [100, 100, 800, 600]);
            hold on;
            
            % Plot band edges for each region
            % ETL
            plot(x(idx_ETL)*1e4, Ec_eV(idx_ETL), 'LineWidth', 2, 'Color', obj.colors{1});
            plot(x(idx_ETL)*1e4, Ev_eV(idx_ETL), 'LineWidth', 2, 'Color', obj.colors{1});
            
            % Absorber
            plot(x(idx_abs)*1e4, Ec_eV(idx_abs), 'LineWidth', 2, 'Color', obj.colors{2});
            plot(x(idx_abs)*1e4, Ev_eV(idx_abs), 'LineWidth', 2, 'Color', obj.colors{2});
            
            % HTL
            plot(x(idx_HTL)*1e4, Ec_eV(idx_HTL), 'LineWidth', 2, 'Color', obj.colors{3});
            plot(x(idx_HTL)*1e4, Ev_eV(idx_HTL), 'LineWidth', 2, 'Color', obj.colors{3});
            
            % Plot quasi-Fermi levels
            plot(x*1e4, Efn_eV, 'k--', 'LineWidth', 1.5);
            plot(x*1e4, Efp_eV, 'k-.', 'LineWidth', 1.5);
            
            % Add region markers
            x_ETL_end = obj.params.L_ETL*1e4;
            x_abs_end = (obj.params.L_ETL + obj.params.L_absorber)*1e4;
            plot([x_ETL_end, x_ETL_end], [min(Ev_eV)-0.5, max(Ec_eV)+0.5], 'k--');
            plot([x_abs_end, x_abs_end], [min(Ev_eV)-0.5, max(Ec_eV)+0.5], 'k--');
            
            % Add labels
            text(x_ETL_end/2, max(Ec_eV)+0.3, 'ETL', 'HorizontalAlignment', 'center');
            text(x_ETL_end + (x_abs_end-x_ETL_end)/2, max(Ec_eV)+0.3, 'Absorber', 'HorizontalAlignment', 'center');
            text(x_abs_end + (obj.params.L_HTL*1e4)/2, max(Ec_eV)+0.3, 'HTL', 'HorizontalAlignment', 'center');
            
            hold off;
            title('Energy Band Diagram');
            xlabel('Position (µm)');
            ylabel('Energy (eV)');
            legend('Conduction Band', 'Valence Band', 'E_{F,n}', 'E_{F,p}', 'Location', 'best');
            grid on;
            
            % Set y-axis limits with some margin
            ylim([min(Ev_eV)-0.5, max(Ec_eV)+0.5]);
        end
        
        function plotCarrierDensities(obj, results)
            % Plot carrier densities with focus on interfaces
            
            % Extract data
            x = obj.params.x;
            n = results.n(end,:)';
            p = results.p(end,:)';
            
            % Create figure
            figure('Name', 'Carrier Densities at Interfaces', 'Position', [100, 100, 1000, 600]);
            
            % Plot carrier densities (log scale)
            semilogy(x*1e4, n, 'b-', 'LineWidth', 2);
            hold on;
            semilogy(x*1e4, p, 'r-', 'LineWidth', 2);
            
            % Highlight interfaces
            x_ETL_end = obj.params.L_ETL*1e4;
            x_abs_end = (obj.params.L_ETL + obj.params.L_absorber)*1e4;
            plot([x_ETL_end, x_ETL_end], [1e0, 1e20], 'k--');
            plot([x_abs_end, x_abs_end], [1e0, 1e20], 'k--');
            
            % Add labels
            text(x_ETL_end/2, 1e19, 'ETL', 'HorizontalAlignment', 'center');
            text(x_ETL_end + (x_abs_end-x_ETL_end)/2, 1e19, 'Absorber', 'HorizontalAlignment', 'center');
            text(x_abs_end + (obj.params.L_HTL*1e4)/2, 1e19, 'HTL', 'HorizontalAlignment', 'center');
            
            hold off;
            title('Carrier Densities Across Interfaces');
            xlabel('Position (µm)');
            ylabel('Carrier Density (cm^{-3})');
            legend('Electrons', 'Holes', 'Location', 'best');
            grid on;
            ylim([1e2, 1e18]);
        end
        
        function plotElectricField(obj, results)
            % Plot electric field across device
            
            % Extract data
            x = obj.params.x;
            E_field = results.E_field;
            
            % Get indices for different regions
            idx_ETL = obj.params.idx_ETL;
            idx_abs = obj.params.idx_absorber;
            idx_HTL = obj.params.idx_HTL;
            
            % Create figure
            figure('Name', 'Electric Field', 'Position', [100, 100, 800, 600]);
            
            % Plot electric field for each region
            hold on;
            plot(x(idx_ETL)*1e4, E_field(idx_ETL)/1e3, 'LineWidth', 2, 'Color', obj.colors{1});
            plot(x(idx_abs)*1e4, E_field(idx_abs)/1e3, 'LineWidth', 2, 'Color', obj.colors{2});
            plot(x(idx_HTL)*1e4, E_field(idx_HTL)/1e3, 'LineWidth', 2, 'Color', obj.colors{3});
            
            % Highlight interfaces
            x_ETL_end = obj.params.L_ETL*1e4;
            x_abs_end = (obj.params.L_ETL + obj.params.L_absorber)*1e4;
            ymin = min(E_field)/1e3 - 5;
            ymax = max(E_field)/1e3 + 5;
            plot([x_ETL_end, x_ETL_end], [ymin, ymax], 'k--');
            plot([x_abs_end, x_abs_end], [ymin, ymax], 'k--');
            
            hold off;
            title('Electric Field');
            xlabel('Position (µm)');
            ylabel('Electric Field (kV/cm)');
            legend('ETL', 'Absorber', 'HTL', 'Location', 'best');
            grid on;
            ylim([ymin, ymax]);
        end
        
        function plotCurrentDensities(obj, results)
            % Plot current densities in device
            
            % Extract data
            x = obj.params.x;
            Jn = results.Jn;
            Jp = results.Jp;
            J_total = results.J_total;
            
            % Get indices for different regions
            idx_ETL = obj.params.idx_ETL;
            idx_abs = obj.params.idx_absorber;
            idx_HTL = obj.params.idx_HTL;
            
            % Create figure
            figure('Name', 'Current Densities', 'Position', [100, 100, 800, 600]);
            
            % Plot current densities
            subplot(2,1,1);
            hold on;
            plot(x*1e4, Jn*1e3, 'b-', 'LineWidth', 2);
            plot(x*1e4, Jp*1e3, 'r-', 'LineWidth', 2);
            plot(x*1e4, J_total*1e3, 'k-', 'LineWidth', 2);
            
            % Highlight interfaces
            x_ETL_end = obj.params.L_ETL*1e4;
            x_abs_end = (obj.params.L_ETL + obj.params.L_absorber)*1e4;
            ymin = min([min(Jn), min(Jp), min(J_total)])*1e3 - 1;
            ymax = max([max(Jn), max(Jp), max(J_total)])*1e3 + 1;
            plot([x_ETL_end, x_ETL_end], [ymin, ymax], 'k--');
            plot([x_abs_end, x_abs_end], [ymin, ymax], 'k--');
            
            hold off;
            title('Current Densities');
            xlabel('Position (µm)');
            ylabel('Current Density (mA/cm^2)');
            legend('Electron', 'Hole', 'Total', 'Location', 'best');
            grid on;
            
            % Plot current components in each region
            subplot(2,1,2);
            
            % Calculate drift and diffusion components
            % For electrons
            Jn_drift = zeros(size(x));
            Jn_diff = zeros(size(x));
            
            for i = 2:length(x)-1
                dx_minus = x(i) - x(i-1);
                dx_plus = x(i+1) - x(i);
                
                % Electron density gradient
                dn_dx = (results.n(end,i+1) - results.n(end,i-1))/(dx_minus + dx_plus);
                
                % Electric field
                E_field = results.E_field(i);
                
                % Calculate components
                Jn_diff(i) = obj.params.q * obj.params.D_n(i) * dn_dx;
                Jn_drift(i) = -obj.params.q * obj.params.mu_n(i) * results.n(end,i) * E_field;
            end
            
            % Plot components
            hold on;
            plot(x*1e4, Jn_drift*1e3, 'b-', 'LineWidth', 2);
            plot(x*1e4, Jn_diff*1e3, 'r-', 'LineWidth', 2);
            plot(x*1e4, (Jn_drift + Jn_diff)*1e3, 'k--', 'LineWidth', 1.5);
            
            % Highlight interfaces
            plot([x_ETL_end, x_ETL_end], [ymin, ymax], 'k--');
            plot([x_abs_end, x_abs_end], [ymin, ymax], 'k--');
            
            hold off;
            title('Electron Current Components');
            xlabel('Position (µm)');
            ylabel('Current Density (mA/cm^2)');
            legend('Drift', 'Diffusion', 'Total', 'Location', 'best');
            grid on;
        end
        
        function plotComparisonWithLight(obj, dark_results, light_results)
            % Compare results with and without illumination
            
            % Extract data
            x = obj.params.x;
            
            % Get indices for different regions
            idx_ETL = obj.params.idx_ETL;
            idx_abs = obj.params.idx_absorber;
            idx_HTL = obj.params.idx_HTL;
            
            % Create figure
            figure('Name', 'Dark vs. Light Comparison', 'Position', [100, 100, 1200, 800]);
            
            % Plot band diagrams
            subplot(2,2,1);
            hold on;
            
            % Convert to eV
            Ec_dark = dark_results.Ec / obj.params.q;
            Ev_dark = dark_results.Ev / obj.params.q;
            Ec_light = light_results.Ec / obj.params.q;
            Ev_light = light_results.Ev / obj.params.q;
            
            % Plot dark conditions
            plot(x*1e4, Ec_dark, 'b-', 'LineWidth', 2);
            plot(x*1e4, Ev_dark, 'r-', 'LineWidth', 2);
            
            % Plot light conditions
            plot(x*1e4, Ec_light, 'b--', 'LineWidth', 1.5);
            plot(x*1e4, Ev_light, 'r--', 'LineWidth', 1.5);
            
            hold off;
            title('Band Diagram: Dark vs. Light');
            xlabel('Position (µm)');
            ylabel('Energy (eV)');
            legend('E_c (dark)', 'E_v (dark)', 'E_c (light)', 'E_v (light)', 'Location', 'best');
            grid on;
            
            % Plot carrier densities
            subplot(2,2,2);
            hold on;
            
            % Dark conditions
            semilogy(x*1e4, dark_results.n(end,:), 'b-', 'LineWidth', 2);
            semilogy(x*1e4, dark_results.p(end,:), 'r-', 'LineWidth', 2);
            
            % Light conditions
            semilogy(x*1e4, light_results.n(end,:), 'b--', 'LineWidth', 1.5);
            semilogy(x*1e4, light_results.p(end,:), 'r--', 'LineWidth', 1.5);
            
            hold off;
            title('Carrier Densities: Dark vs. Light');
            xlabel('Position (µm)');
            ylabel('Carrier Density (cm^{-3})');
            legend('n (dark)', 'p (dark)', 'n (light)', 'p (light)', 'Location', 'best');
            grid on;
            ylim([1e2, 1e18]);
            
            % Plot quasi-Fermi levels
            subplot(2,2,3);
            hold on;
            
            % Convert to eV
            Efn_dark = dark_results.Efn / obj.params.q;
            Efp_dark = dark_results.Efp / obj.params.q;
            Efn_light = light_results.Efn / obj.params.q;
            Efp_light = light_results.Efp / obj.params.q;
            
            % Plot quasi-Fermi levels
            plot(x*1e4, Efn_dark, 'b-', 'LineWidth', 2);
            plot(x*1e4, Efp_dark, 'r-', 'LineWidth', 2);
            plot(x*1e4, Efn_light, 'b--', 'LineWidth', 1.5);
            plot(x*1e4, Efp_light, 'r--', 'LineWidth', 1.5);
            
            hold off;
            title('Quasi-Fermi Levels: Dark vs. Light');
            xlabel('Position (µm)');
            ylabel('Energy (eV)');
            legend('E_{F,n} (dark)', 'E_{F,p} (dark)', 'E_{F,n} (light)', 'E_{F,p} (light)', 'Location', 'best');
            grid on;
            
            % Plot current densities
            subplot(2,2,4);
            hold on;
            
            % Dark conditions
            plot(x*1e4, dark_results.J_total*1e3, 'k-', 'LineWidth', 2);
            
            % Light conditions
            plot(x*1e4, light_results.J_total*1e3, 'r-', 'LineWidth', 2);
            
            hold off;
            title('Total Current: Dark vs. Light');
            xlabel('Position (µm)');
            ylabel('Current Density (mA/cm^2)');
            legend('Dark', 'Light', 'Location', 'best');
            grid on;
        end
    end
end