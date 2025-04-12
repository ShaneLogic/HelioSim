classdef InterfaceHandler < handle
    % InterfaceHandler - Manages material interfaces in solar cell
    % Ensures proper continuity conditions across heterojunctions
    
    properties
        params      % Solar cell parameters
    end
    
    methods
        function obj = InterfaceHandler(params)
            % Constructor - initialize with cell parameters
            obj.params = params;
        end
        
        function [n_interface, p_interface] = applyConditions(obj, phi, n, p, idx)
            % Apply interface conditions at a specified index
            
            % Determine which interface we're at
            if idx == obj.params.idx_ETL(end)
                % ETL/Absorber interface
                interface_type = 'ETL_abs';
            else
                % Absorber/HTL interface
                interface_type = 'abs_HTL';
            end
            
            % Get material parameters on both sides
            if strcmp(interface_type, 'ETL_abs')
                % ETL side (left)
                chi_L = obj.params.chi_ETL;
                Eg_L = obj.params.Eg_ETL;
                Nc_L = obj.params.Nc_ETL;
                Nv_L = obj.params.Nv_ETL;
                
                % Absorber side (right)
                chi_R = obj.params.chi_abs;
                Eg_R = obj.params.Eg_abs;
                Nc_R = obj.params.Nc_abs;
                Nv_R = obj.params.Nv_abs;
            else
                % Absorber side (left)
                chi_L = obj.params.chi_abs;
                Eg_L = obj.params.Eg_abs;
                Nc_L = obj.params.Nc_abs;
                Nv_L = obj.params.Nv_abs;
                
                % HTL side (right)
                chi_R = obj.params.chi_HTL;
                Eg_R = obj.params.Eg_HTL;
                Nc_R = obj.params.Nc_HTL;
                Nv_R = obj.params.Nv_HTL;
            end
            
            % Calculate band offsets
            q = obj.params.q;
            deltaEc = q * (chi_L - chi_R);
            deltaEv = deltaEc - q * (Eg_L - Eg_R);
            
            % Get thermal energy
            kT = obj.params.kb * obj.params.T;
            
            % Get neighboring indices
            idx_L = idx - 1;
            idx_R = idx + 1;
            
            % Get potential values
            phi_L = phi(idx_L);
            phi_R = phi(idx_R);
            phi_I = phi(idx);
            
            % Calculate band energies at the interface and adjacent points
            Ec_L = -q * chi_L - q * phi_L;
            Ec_R = -q * chi_R - q * phi_R;
            Ec_I = -q * chi_L - q * phi_I; % Reference to left material
            
            Ev_L = Ec_L - q * Eg_L;
            Ev_R = Ec_R - q * Eg_R;
            Ev_I = Ec_I - q * Eg_L; % Reference to left material
            
            % Calculate quasi-Fermi levels on left and right
            Efn_L = Ec_L + kT * log(n(idx_L) / Nc_L);
            Efn_R = Ec_R + kT * log(n(idx_R) / Nc_R);
            
            Efp_L = Ev_L + kT * log(p(idx_L) / Nv_L);
            Efp_R = Ev_R + kT * log(p(idx_R) / Nv_R);
            
            % Enforce continuity of quasi-Fermi levels at interface
            % Average approach for numerical stability
            Efn_I = (Efn_L + Efn_R) / 2;
            Efp_I = (Efp_L + Efp_R) / 2;
            
            % Calculate carrier densities at interface
            n_interface = Nc_L * exp((Efn_I - Ec_I) / kT);
            p_interface = Nv_L * exp((Efp_I - Ev_I) / kT);
                        % Handle special cases for numerical stability
            % Ensure carrier densities don't get too small
            n_min = 1e2;
            p_min = 1e2;
            n_interface = max(n_interface, n_min);
            p_interface = max(p_interface, p_min);
        end
    end
end