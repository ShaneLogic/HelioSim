classdef OpticalGeneration < handle
    % OpticalGeneration - Class to calculate optical generation profiles
    % Models light absorption and carrier generation in the solar cell
    
    properties
        params              % Solar cell parameters
        AM15                % AM1.5 solar spectrum data
    end
    
    methods
        function obj = OpticalGeneration(params)
            % Constructor - initialize with solar cell parameters
            obj.params = params;
            
            % Load or create AM1.5 spectrum data
            obj.loadAM15Spectrum();
        end
        
        function loadAM15Spectrum(obj)
            % Load AM1.5G solar spectrum data
            % If no file available, use approximate model
            
            try
                % Try to load spectrum data from file
                % Format: [wavelength(nm), irradiance(W/m^2/nm)]
                load('AM15G.mat', 'spectrum');
                obj.AM15 = spectrum;
            catch
                % If file not found, create simple approximation
                wavelength = (300:5:1200)';  % nm
                % Simple approximation of AM1.5G
                irradiance = 1000/900 * exp(-(wavelength-700).^2/160000);
                obj.AM15 = [wavelength, irradiance];
            end
        end
        
        function G = calculateGeneration(obj)
            % Calculate generation rate profile across the device
            
            % Get position grid
            x = obj.params.x;
            Nx = length(x);
            
            % Initialize generation rate array
            G = zeros(Nx, 1);
            
            % Only generate carriers in the absorber layer
            idx_abs = obj.params.idx_absorber;
            
            % Calculate generation profile in absorber
            % Simple exponential decay model (Beer-Lambert law)
            x_abs = x(idx_abs);
            x_front = x_abs(1);  % Front of absorber
            
            % Account for front surface reflection
            transmitted = 1 - obj.params.R_front;
            
            % Apply Beer-Lambert law: G(x) = G₀ * (1-R) * exp(-α(x-x₀))
            G(idx_abs) = obj.params.G_max * transmitted * ...
                          exp(-obj.params.alpha_abs * (x_abs - x_front));
            
            % Add some small generation in ETL (UV absorption)
            idx_ETL = obj.params.idx_ETL;
            G(idx_ETL) = 0.05 * obj.params.G_max * transmitted * ...
                          exp(-5e5 * (x(idx_ETL) - x(1)));
            
            % Very little generation in HTL
            idx_HTL = obj.params.idx_HTL;
            G(idx_HTL) = 0.001 * obj.params.G_max * transmitted * ...
                          exp(-obj.params.alpha_abs * (x(idx_HTL) - x_front));
            
            % If device is not illuminated, return zeros
            if ~obj.params.illumination
                G = zeros(Nx, 1);
            end
        end
        
        function G_detailed = calculateDetailedGeneration(obj)
            % Calculate wavelength-resolved generation for detailed analysis
            % Uses full AM1.5 spectrum and wavelength-dependent absorption
            
            % Get position grid
            x = obj.params.x;
            Nx = length(x);
            
            % Get wavelength grid from spectrum
            wavelengths = obj.AM15(:,1);  % nm
            irradiance = obj.AM15(:,2);   % W/m^2/nm
            
            % Initialize array for generation rate
            G_detailed = zeros(Nx, length(wavelengths));
            
            % Only consider absorber layer
            idx_abs = obj.params.idx_absorber;
            x_abs = x(idx_abs);
            x_front = x_abs(1);
            
            % Calculate wavelength-dependent absorption coefficient
            alpha_wavelength = obj.calculateAlpha(wavelengths);
            
            % Account for reflection
            transmitted = 1 - obj.params.R_front;
            
            % Calculate generation at each position and wavelength
            for i = 1:length(wavelengths)
                lambda = wavelengths(i);
                alpha = alpha_wavelength(i);
                
                % Photon energy (J)
                E_photon = 1.24e-6 / (lambda * 1e-9);
                
                % Photon flux (#/m^2/s/nm)
                photon_flux = irradiance(i) / E_photon;
                
                % Generation rate in absorber (#/cm^3/s)
                G_detailed(idx_abs, i) = transmitted * photon_flux * alpha * ...
                                         exp(-alpha * (x_abs - x_front)) * 1e-4;
            end
            
            % Total generation rate at each position
            G_total = sum(G_detailed, 2);
            
            % Normalize to match maximum generation rate
            G_total = G_total * (obj.params.G_max / max(G_total));
            
            % If not illuminated, return zeros
            if ~obj.params.illumination
                G_detailed = zeros(size(G_detailed));
            end
        end
        
        function alpha = calculateAlpha(obj, wavelengths)
            % Calculate wavelength-dependent absorption coefficient
            % Simplified model based on empirical data
            
            % Initialize alpha
            alpha = zeros(size(wavelengths));
            
            % Calculate absorption edge wavelength from bandgap
            lambda_g = 1240 / obj.params.Eg_abs;  % nm
            
            % Simple model:
            % 1. Below bandgap: weak absorption
            % 2. Near bandgap: sqrt dependence
            % 3. Above bandgap: strong absorption
            for i = 1:length(wavelengths)
                lambda = wavelengths(i);
                if lambda >= lambda_g
                    % Below bandgap - weak absorption
                    alpha(i) = obj.params.alpha_abs * 0.01 * exp(-(lambda-lambda_g)/50);
                else
                    % Above bandgap - strong absorption
                    % Square root dependency near band edge
                    alpha(i) = obj.params.alpha_abs * sqrt(lambda_g/lambda - 1) + 0.1*obj.params.alpha_abs;
                end
            end
        end
    end
end