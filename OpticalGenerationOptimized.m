classdef OpticalGenerationOptimized < handle
    % OpticalGenerationOptimized - 优化的光生载流子生成模型
    % 实现波长依赖的光吸收和载流子生成，针对钙钛矿材料进行了优化
    
    properties
        params              % 太阳能电池参数
        AM15                % AM1.5太阳光谱数据
        enable_interference % 是否启用干涉效应
        enable_back_reflection % 是否启用背面反射
    end
    
    methods
        function obj = OpticalGenerationOptimized(params)
            % 构造函数 - 初始化光生载流子生成模型
            obj.params = params;
            
            % 加载AM1.5光谱数据
            obj.loadAM15Spectrum();
            
            % 默认不启用干涉效应和背面反射
            obj.enable_interference = false;
            obj.enable_back_reflection = true;
        end
        
        function loadAM15Spectrum(obj)
            % 加载AM1.5G太阳光谱数据
            % 如果没有文件，使用近似模型
            
            try
                % 尝试从文件加载光谱数据
                % 格式: [wavelength(nm), irradiance(W/m^2/nm)]
                load('AM15G.mat', 'spectrum');
                obj.AM15 = spectrum;
            catch
                % 如果文件未找到，创建简单近似
                wavelength = (300:5:1200)';  % nm
                % AM1.5G的简单近似
                irradiance = 1000/900 * exp(-(wavelength-700).^2/160000);
                obj.AM15 = [wavelength, irradiance];
            end
        end
        
        function G = calculateGeneration(obj)
            % 计算光生载流子生成率分布
            % 使用简化的Beer-Lambert吸收模型
            
            % 获取位置网格
            x = obj.params.x;
            Nx = length(x);
            
            % 初始化生成率数组
            G = zeros(Nx, 1);
            
            % 如果没有光照，返回零生成率
            if ~obj.params.illumination
                return;
            end
            
            % 仅在吸收层生成载流子
            idx_abs = obj.params.idx_absorber;
            
            % 计算吸收层中的生成率分布
            % 使用Beer-Lambert定律: G(x) = G₀ * (1-R) * exp(-α(x-x₀))
            x_abs = x(idx_abs);
            x_front = x_abs(1);  % 吸收层前表面
            
            % 考虑前表面反射
            transmitted = 1 - obj.params.R_front;
            
            if obj.enable_interference
                % 使用干涉模型计算生成率
                G(idx_abs) = obj.calculateInterferenceGeneration(x_abs);
            else
                % 使用改进的Beer-Lambert定律计算光生率
                % 使用带隔的吸收系数，更符合MAPbI3的实际吸收特性
                % 吸收系数与波长相关，这里使用平均值
                
                % 计算光强度分布 - 前向传播
                I_forward = obj.params.G_max * transmitted * ...
                            exp(-obj.params.alpha_abs * (x_abs - x_front));
                
                % 计算光生率 - 前向传播部分
                G_forward = obj.params.alpha_abs * I_forward;
                
                % 如果启用背面反射，考虑反向传播的光
                if obj.enable_back_reflection
                    % 使用更准确的背面反射率
                    R_back = 0.8;  % 背面反射率，考虑到HTL/吸收层界面的折射率
                    x_back = x_abs(end);  % 吸收层后表面
                    
                    % 计算到达后表面的光强度
                    I_back = obj.params.G_max * transmitted * exp(-obj.params.alpha_abs * (x_back - x_front));
                    
                    % 计算反射后的光强度分布 - 反向传播
                    I_backward = I_back * R_back * exp(-obj.params.alpha_abs * (x_back - x_abs));
                    
                    % 计算光生率 - 反向传播部分
                    G_backward = obj.params.alpha_abs * I_backward;
                    
                    % 总光生率 = 前向 + 反向
                    G(idx_abs) = G_forward + G_backward;
                else
                    % 如果不考虑背面反射，只使用前向传播部分
                    G(idx_abs) = G_forward;
                end
            end
            
            % 在ETL中添加光生生成 (UV吸收)
            idx_ETL = obj.params.idx_ETL;
            
            % ETL的吸收系数（TiO2对UV光有较强吸收）
            alpha_ETL = 5e4; % [cm^-1]
            
            % 计算ETL中的光强度分布
            I_ETL = obj.params.G_max * transmitted * exp(-alpha_ETL * (x(idx_ETL) - x(1)));
            
            % 计算ETL中的光生率
            G(idx_ETL) = 0.1 * alpha_ETL * I_ETL; % 系数调整以反映TiO2对可见光的较低吸收
            
            % HTL中的光生生成
            idx_HTL = obj.params.idx_HTL;
            
            % 如果启用背面反射，考虑到达吸收层/HTL界面的光强度
            if obj.enable_back_reflection
                % 计算到达吸收层后表面的光强度
                x_back = x(obj.params.idx_absorber(end));
                I_back = obj.params.G_max * transmitted * exp(-obj.params.alpha_abs * (x_back - x_front));
                
                % HTL的吸收系数（Spiro-OMeTAD对可见光有较低吸收）
                alpha_HTL = 1e3; % [cm^-1]
                
                % 计算界面反射系数
                n_abs = 2.5 + 0.5i;  % MAPbI3的复折射率
                n_HTL = 1.8;  % Spiro-OMeTAD的折射率
                r_abs_HTL = (n_abs - n_HTL)/(n_abs + n_HTL);  % 吸收层/HTL界面反射系数
                R_abs_HTL = abs(r_abs_HTL)^2;  % 界面反射率
                
                % 计算HTL中的光强度分布 - 使用计算的反射系数
                I_HTL = I_back * (1 - R_abs_HTL) * exp(-alpha_HTL * (x(idx_HTL) - x_back));
                
                % 计算HTL中的光生率
                G(idx_HTL) = 0.05 * alpha_HTL * I_HTL; % 系数调整以反映Spiro-OMeTAD的量子效率
            else
                % 如果不考虑背面反射，则HTL中几乎没有光生生成
                G(idx_HTL) = 0.001 * obj.params.G_max * transmitted * exp(-obj.params.alpha_abs * (x(idx_HTL) - x_front));
            end
        end
        
        function G_interference = calculateInterferenceGeneration(obj, x_abs)
            % 计算考虑干涉效应的生成率
            % 使用传输矩阵法计算多层结构中的光场分布
            
            % 获取波长网格
            wavelengths = obj.AM15(:,1);  % nm
            irradiance = obj.AM15(:,2);   % W/m^2/nm
            
            % 初始化生成率数组
            G_interference = zeros(length(x_abs), 1);
            
            % 计算每个波长的贡献
            for i = 1:10:length(wavelengths)  % 每10个波长点计算一次，提高效率
                lambda = wavelengths(i);
                
                % 计算复折射率
                n_ETL = 2.5;  % TiO2的折射率
                n_abs = 2.5 + 0.5i;  % MAPbI3的复折射率
                n_HTL = 1.8;  % Spiro-OMeTAD的折射率
                
                % 计算界面反射系数
                r_ETL_abs = (n_ETL - n_abs)/(n_ETL + n_abs);  % ETL/吸收层界面反射系数
                r_abs_HTL = (n_abs - n_HTL)/(n_abs + n_HTL);  % 吸收层/HTL界面反射系数
                
                % 计算波数
                k0 = 2*pi/(lambda*1e-9);
                
                % 计算光场分布
                E_field = ones(length(x_abs), 1);
                
                % 简化的干涉模型：正向和反向传播波的叠加
                for j = 1:length(x_abs)
                    x_rel = x_abs(j) - x_abs(1);  % 相对于吸收层前表面的位置
                    
                    % 正向传播波
                    E_forward = exp(-1i * k0 * n_abs * x_rel);
                    
                    % 反向传播波 (来自背面反射)
                    % 使用计算的界面反射系数而不是硬编码的值
                    R_back = abs(r_abs_HTL)^2;  % 使用界面反射系数计算反射率
                    x_total = x_abs(end) - x_abs(1);  % 吸收层总厚度
                    E_backward = R_back * exp(-1i * k0 * n_abs * (2*x_total - x_rel));
                    
                    % 总场
                    E_field(j) = abs(E_forward + E_backward)^2;
                end
                
                % 计算吸收
                alpha = 4*pi*imag(n_abs)/(lambda*1e-9);  % 吸收系数
                absorption = alpha * E_field;
                
                % 计算光子能量 (J)
                E_photon = obj.params.h * obj.params.c / (lambda * 1e-9);
                
                % 计算光子通量 (#/m^2/s/nm)
                photon_flux = irradiance(i) / E_photon;
                
                % 生成率贡献 (#/cm^3/s)
                G_contribution = absorption .* photon_flux * 1e-4;
                
                % 添加到总生成率
                G_interference = G_interference + G_contribution;
            end
            
            % 归一化到最大生成率
            G_interference = G_interference * (obj.params.G_max / max(G_interference));
        end
        
        function G_detailed = calculateDetailedGeneration(obj)
            % 计算波长分辨的生成率，用于详细分析
            % 使用完整的AM1.5光谱和波长依赖的吸收
            
            % 获取位置网格
            x = obj.params.x;
            Nx = length(x);
            
            % 获取波长网格
            wavelengths = obj.AM15(:,1);  % nm
            irradiance = obj.AM15(:,2);   % W/m^2/nm
            
            % 初始化生成率数组
            G_detailed = zeros(Nx, length(wavelengths));
            
            % 如果没有光照，返回零生成率
            if ~obj.params.illumination
                return;
            end
            
            % 仅考虑吸收层
            idx_abs = obj.params.idx_absorber;
            x_abs = x(idx_abs);
            x_front = x_abs(1);
            
            % 计算波长依赖的吸收系数
            alpha_wavelength = obj.calculateAlpha(wavelengths);
            
            % 考虑反射
            transmitted = 1 - obj.params.R_front;
            
            % 计算每个位置和波长的生成率
            for i = 1:length(wavelengths)
                lambda = wavelengths(i);
                alpha = alpha_wavelength(i);
                
                % 光子能量 (J)
                E_photon = obj.params.h * obj.params.c / (lambda * 1e-9);
                
                % 光子通量 (#/m^2/s/nm)
                photon_flux = irradiance(i) / E_photon;
                
                % 吸收层中的生成率 (#/cm^3/s)
                G_detailed(idx_abs, i) = transmitted * photon_flux * alpha * ...
                                         exp(-alpha * (x_abs - x_front)) * 1e-4;
                
                % 如果启用背面反射，添加反向传播的光
                if obj.enable_back_reflection
                    R_back = 0.9;  % 背面反射率
                    x_back = x_abs(end);  % 吸收层后表面
                    G_reflected = transmitted * photon_flux * alpha * R_back * ...
                                  exp(-alpha * (2*x_back - x_abs - x_front)) * 1e-4;
                    G_detailed(idx_abs, i) = G_detailed(idx_abs, i) + G_reflected;
                end
            end
            
            % 如果没有光照，返回零生成率
            if ~obj.params.illumination
                G_detailed = zeros(size(G_detailed));
            end
        end
        
        function alpha = calculateAlpha(obj, wavelengths)
            % 计算波长依赖的吸收系数
            % 基于钙钛矿材料的经验模型
            
            % 初始化alpha
            alpha = zeros(size(wavelengths));
            
            % 从带隙计算吸收边缘波长
            lambda_g = 1240 / obj.params.Eg_abs;  % nm
            
            % 简单模型:
            % 1. 带隙以下: 弱吸收
            % 2. 带隙附近: 平方根依赖
            % 3. 带隙以上: 强吸收
            for i = 1:length(wavelengths)
                lambda = wavelengths(i);
                if lambda >= lambda_g
                    % 带隙以下 - 弱吸收
                    alpha(i) = obj.params.alpha_abs * 0.01 * exp(-(lambda-lambda_g)/50);
                else
                    % 带隙以上 - 强吸收
                    % 带边附近的平方根依赖
                    alpha(i) = obj.params.alpha_abs * sqrt(lambda_g/lambda - 1) + 0.1*obj.params.alpha_abs;
                end
            end
        end
        
        function setInterferenceModel(obj, enable)
            % 设置是否启用干涉模型
            obj.enable_interference = enable;
        end
        
        function setBackReflection(obj, enable)
            % 设置是否启用背面反射
            obj.enable_back_reflection = enable;
        end
        
        function plotAbsorptionSpectrum(obj)
            % 绘制吸收光谱
            
            % 计算波长依赖的吸收系数
            wavelengths = obj.AM15(:,1);  % nm
            alpha = obj.calculateAlpha(wavelengths);
            
            % 计算吸收率 (1-exp(-alpha*d))
            d = obj.params.L_absorber;  % 吸收层厚度
            absorption = 1 - exp(-alpha * d);
            
            % 绘制吸收光谱
            figure('Name', '吸收光谱');
            
            % 吸收系数
            subplot(2,1,1);
            semilogy(wavelengths, alpha, 'LineWidth', 2);
            xlabel('波长 (nm)');
            ylabel('吸收系数 (cm^{-1})');
            title('钙钛矿吸收系数');
            grid on;
            xlim([300, 1000]);
            
            % 吸收率
            subplot(2,1,2);
            plot(wavelengths, absorption * 100, 'LineWidth', 2);
            xlabel('波长 (nm)');
            ylabel('吸收率 (%)');
            title(['吸收层厚度: ' num2str(d*1e7) ' nm']);
            grid on;
            xlim([300, 1000]);
            ylim([0, 100]);
        end
    end
end
