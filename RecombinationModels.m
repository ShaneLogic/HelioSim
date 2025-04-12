classdef RecombinationModels < handle
    % RecombinationModels - Class to calculate recombination rates
    % Includes SRH, radiative, and Auger recombination mechanisms
    
    properties
        params                  % Solar cell parameters
        
        % Recombination parameters
        B_rad_abs = 1e-10;     % Radiative recombination coefficient [cm^3/s]
        Cn_auger_abs = 1e-30;  % Auger coefficient (electrons) [cm^6/s]
        Cp_auger_abs = 1e-30;  % Auger coefficient (holes) [cm^6/s]
    end
    
    methods
        function obj = RecombinationModels(params)
            % Constructor - initialize with solar cell parameters
            obj.params = params;
        end
        
        function R = calculateRecombination(obj, n, p)
            % Calculate recombination rates at each position
            
            % Get total number of points
            Nx = obj.params.Nx_total;
            
            % Initialize recombination rate array
            R = zeros(Nx, 1);
            
            % Calculate recombination in each region
            % ETL region - only SRH
            idx = obj.params.idx_ETL;
            R(idx) = obj.SRH_recombination(n(idx), p(idx), ...
                      obj.params.ni_ETL, obj.params.tau_n_ETL, obj.params.tau_p_ETL);
            
            % Absorber region - SRH, radiative, and Auger
            idx = obj.params.idx_absorber;
            R_SRH = obj.SRH_recombination(n(idx), p(idx), ...
                    obj.params.ni_abs, obj.params.tau_n_abs, obj.params.tau_p_abs);
            R_rad = obj.radiative_recombination(n(idx), p(idx), obj.params.ni_abs);
            R_auger = obj.Auger_recombination(n(idx), p(idx), obj.params.ni_abs);
            R(idx) = R_SRH + R_rad + R_auger;
            
            % HTL region - only SRH
            idx = obj.params.idx_HTL;
            R(idx) = obj.SRH_recombination(n(idx), p(idx), ...
                      obj.params.ni_HTL, obj.params.tau_n_HTL, obj.params.tau_p_HTL);
        end
        
        function R = SRH_recombination(obj, n, p, ni, tau_n, tau_p)
            % Calculate Shockley-Read-Hall recombination rate
            % R_SRH = (np - ni²) / (τ_n(p + ni) + τ_p(n + ni))
            
            % Calculate recombination rate
            R = (n .* p - ni^2) ./ (tau_n * (p + ni) + tau_p * (n + ni));
            
            % Handle numerical issues
            R(isnan(R)) = 0;
            R(isinf(R)) = 0;
        end
        
        function R = radiative_recombination(obj, n, p, ni)
            % Calculate radiative (band-to-band) recombination rate
            % R_rad = B(np - ni²)
            
            % Calculate recombination rate
            R = obj.B_rad_abs * (n .* p - ni^2);
            
            % Handle numerical issues
            R(isnan(R)) = 0;
            R(isinf(R)) = 0;
        end
        
        function R = Auger_recombination(obj, n, p, ni)
            % Calculate Auger recombination rate
            % R_Auger = (np - ni²)(C_n n + C_p p)
            
            % Calculate recombination rate
            R = (n .* p - ni^2) .* (obj.Cn_auger_abs * n + obj.Cp_auger_abs * p);
            
            % Handle numerical issues
            R(isnan(R)) = 0;
            R(isinf(R)) = 0;
        end
        
        function R_total = calculateTotalRecombination(obj, results)
            % Calculate total recombination rate for post-processing
            
            % Get steady-state carrier densities
            n = results.n(end,:)';
            p = results.p(end,:)';
            
            % Calculate recombination
            R = obj.calculateRecombination(n, p);
            
            % Calculate total volumetric recombination
            dx = diff(obj.params.x);
            dx_full = [dx(1); (dx(1:end-1) + dx(2:end))/2; dx(end)];
            
            % Integrate recombination rate over device
            R_total = sum(R .* dx_full);
        end
    end
end