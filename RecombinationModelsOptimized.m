classdef RecombinationModelsOptimized < handle
    % RecombinationModelsOptimized - 优化的复合模型
    % 实现SRH、Auger和辐射复合，针对钙钛矿太阳能电池进行了优化
    
    properties
        params          % 太阳能电池参数
        
        % SRH复合参数
        tau_n           % 电子SRH寿命 [s]
        tau_p           % 空穴SRH寿命 [s]
        Et              % 陷阱能级 (相对于本征费米能级) [eV]
        
        % Auger复合参数
        Cn              % 电子Auger系数 [cm^6/s]
        Cp              % 空穴Auger系数 [cm^6/s]
        
        % 辐射复合参数
        B               % 辐射复合系数 [cm^3/s]
        
        % 界面复合参数
        S_ETL_abs       % ETL/吸收层界面复合速率 [cm/s]
        S_abs_HTL       % 吸收层/HTL界面复合速率 [cm/s]
        
        % 位置依赖的参数
        position_dependent % 是否使用位置依赖的参数
    end
    
    methods
        function obj = RecombinationModelsOptimized(params)
            % 构造函数 - 初始化复合模型
            obj.params = params;
            
            % 计算本征载流子浓度
            % 对于每个材料，ni = sqrt(Nc*Nv)*exp(-Eg/(2*kT))
            % ETL层
            params.ni_ETL = sqrt(params.Nc_ETL * params.Nv_ETL) * exp(-params.Eg_ETL * params.q / (2 * params.kb * params.T));
            % 吸收层
            params.ni_abs = sqrt(params.Nc_abs * params.Nv_abs) * exp(-params.Eg_abs * params.q / (2 * params.kb * params.T));
            % HTL层
            params.ni_HTL = sqrt(params.Nc_HTL * params.Nv_HTL) * exp(-params.Eg_HTL * params.q / (2 * params.kb * params.T));
            
            % 设置默认参数
            obj.setDefaultParameters();
            
            % 默认不使用位置依赖的参数
            obj.position_dependent = false;
        end
        
        function setDefaultParameters(obj)
            % 设置默认复合参数
            
            % SRH复合
            obj.tau_n = 1e-6;  % 1 μs
            obj.tau_p = 1e-6;  % 1 μs
            obj.Et = 0;        % 中间带隙
            
            % Auger复合
            obj.Cn = 1e-30;    % 较小的Auger系数
            obj.Cp = 1e-30;
            
            % 辐射复合
            obj.B = 1e-10;     % 典型的辐射复合系数
            
            % 界面复合
            obj.S_ETL_abs = 1e4;  % 10^4 cm/s
            obj.S_abs_HTL = 1e4;  % 10^4 cm/s
        end
        
        function setPerovskiteParameters(obj)
            % 设置钙钛矿太阳能电池的复合参数
            
            % SRH复合 (钙钛矿通常有较长的载流子寿命)
            obj.tau_n = 1e-6;  % 1 μs
            obj.tau_p = 1e-6;  % 1 μs
            obj.Et = 0;        % 中间带隙
            
            % Auger复合 (钙钛矿中通常不显著)
            obj.Cn = 1e-31;
            obj.Cp = 1e-31;
            
            % 辐射复合 (钙钛矿有较强的辐射复合)
            obj.B = 5e-10;
            
            % 界面复合 (钙钛矿器件中的关键损失机制)
            obj.S_ETL_abs = 1e4;  % TiO2/MAPbI3界面
            obj.S_abs_HTL = 1e4;  % MAPbI3/Spiro-OMeTAD界面
        end
        
        function enablePositionDependentParameters(obj, enable)
            % 启用或禁用位置依赖的参数
            obj.position_dependent = enable;
        end
        
        function R = calculateTotalRecombination(obj, n, p, x)
            % 计算总复合率
            % n: 电子密度 [cm^-3]
            % p: 空穴密度 [cm^-3]
            % x: 位置网格 [cm]
            
            % 初始化总复合率
            R = zeros(size(n));
            
            % 计算体复合
            R_SRH = obj.calculateSRHRecombination(n, p, x);
            R_Auger = obj.calculateAugerRecombination(n, p, x);
            R_rad = obj.calculateRadiativeRecombination(n, p, x);
            
            % 总体复合率
            R = R_SRH + R_Auger + R_rad;
            
            % 添加界面复合 (在界面位置的网格点上)
            R = obj.addInterfaceRecombination(R, n, p, x);
        end
        
        function R_SRH = calculateSRHRecombination(obj, n, p, x)
            % 计算SRH复合率
            % R_SRH = (np - ni^2) / (tau_n*(n+n1) + tau_p*(p+p1))
            
            % 获取本征载流子密度
            % 根据位置选择不同材料的本征载流子密度
            ni = zeros(size(x));
            
            % 确定每个点所属的材料层
            for i = 1:length(x)
                if x(i) <= obj.params.L_ETL
                    % ETL层
                    ni(i) = obj.params.ni_ETL;
                elseif x(i) <= obj.params.L_ETL + obj.params.L_absorber
                    % 吸收层
                    ni(i) = obj.params.ni_abs;
                else
                    % HTL层
                    ni(i) = obj.params.ni_HTL;
                end
            end
            
            % 计算n1和p1 (与陷阱能级相关的参数)
            % 根据位置选择不同材料的能带参数
            Ec = zeros(size(x));
            Ev = zeros(size(x));
            
            % 确定每个点所属的材料层
            for i = 1:length(x)
                if x(i) <= obj.params.L_ETL
                    % ETL层
                    Ec(i) = obj.params.Ec_ETL;
                    Ev(i) = obj.params.Ev_ETL;
                elseif x(i) <= obj.params.L_ETL + obj.params.L_absorber
                    % 吸收层
                    Ec(i) = obj.params.Ec_abs;
                    Ev(i) = obj.params.Ev_abs;
                else
                    % HTL层
                    Ec(i) = obj.params.Ec_HTL;
                    Ev(i) = obj.params.Ev_HTL;
                end
            end
            
            % 计算本征能级
            Ei = (Ec + Ev) / 2;  % 本征能级
            E_trap = Ei + obj.params.q * obj.Et;     % 陷阱能级
            
            % 根据位置计算n1和p1
            n1 = zeros(size(x));
            p1 = zeros(size(x));
            
            for i = 1:length(x)
                if x(i) <= obj.params.L_ETL
                    % ETL层
                    n1(i) = obj.params.Nc_ETL * exp(-(Ec(i) - E_trap(i)) / (obj.params.kb * obj.params.T));
                    p1(i) = obj.params.Nv_ETL * exp(-(E_trap(i) - Ev(i)) / (obj.params.kb * obj.params.T));
                elseif x(i) <= obj.params.L_ETL + obj.params.L_absorber
                    % 吸收层
                    n1(i) = obj.params.Nc_abs * exp(-(Ec(i) - E_trap(i)) / (obj.params.kb * obj.params.T));
                    p1(i) = obj.params.Nv_abs * exp(-(E_trap(i) - Ev(i)) / (obj.params.kb * obj.params.T));
                else
                    % HTL层
                    n1(i) = obj.params.Nc_HTL * exp(-(Ec(i) - E_trap(i)) / (obj.params.kb * obj.params.T));
                    p1(i) = obj.params.Nv_HTL * exp(-(E_trap(i) - Ev(i)) / (obj.params.kb * obj.params.T));
                end
            end
            
            % 获取SRH寿命
            tau_n = zeros(size(x));
            tau_p = zeros(size(x));
            
            % 根据位置设置寿命
            for i = 1:length(x)
                if x(i) <= obj.params.L_ETL
                    % ETL层
                    tau_n(i) = obj.params.tau_n_ETL;
                    tau_p(i) = obj.params.tau_p_ETL;
                elseif x(i) <= obj.params.L_ETL + obj.params.L_absorber
                    % 吸收层
                    tau_n(i) = obj.params.tau_n_abs;
                    tau_p(i) = obj.params.tau_p_abs;
                else
                    % HTL层
                    tau_n(i) = obj.params.tau_n_HTL;
                    tau_p(i) = obj.params.tau_p_HTL;
                end
            end
            
            % 计算SRH复合率
            R_SRH = (n .* p - ni.^2) ./ (tau_n .* (n + n1) + tau_p .* (p + p1));
        end
        
        function R_Auger = calculateAugerRecombination(obj, n, p, x)
            % 计算Auger复合率
            % R_Auger = Cn*n*(np - ni^2) + Cp*p*(np - ni^2)
            
            % 获取本征载流子密度
            % 根据位置选择不同材料的本征载流子密度
            ni = zeros(size(x));
            
            % 确定每个点所属的材料层
            for i = 1:length(x)
                if x(i) <= obj.params.L_ETL
                    % ETL层
                    ni(i) = obj.params.ni_ETL;
                elseif x(i) <= obj.params.L_ETL + obj.params.L_absorber
                    % 吸收层
                    ni(i) = obj.params.ni_abs;
                else
                    % HTL层
                    ni(i) = obj.params.ni_HTL;
                end
            end
            
            % 获取Auger系数
            if obj.position_dependent
                [Cn, Cp] = obj.getPositionDependentAugerCoefficients(x);
            else
                Cn = obj.Cn * ones(size(x));
                Cp = obj.Cp * ones(size(x));
            end
            
            % 计算Auger复合率
            R_Auger = Cn .* n .* (n .* p - ni.^2) + Cp .* p .* (n .* p - ni.^2);
        end
        
        function R_rad = calculateRadiativeRecombination(obj, n, p, x)
            % 计算辐射复合率
            % R_rad = B*(np - ni^2)
            
            % 获取本征载流子密度
            % 根据位置选择不同材料的本征载流子密度
            ni = zeros(size(x));
            
            % 确定每个点所属的材料层
            for i = 1:length(x)
                if x(i) <= obj.params.L_ETL
                    % ETL层
                    ni(i) = obj.params.ni_ETL;
                elseif x(i) <= obj.params.L_ETL + obj.params.L_absorber
                    % 吸收层
                    ni(i) = obj.params.ni_abs;
                else
                    % HTL层
                    ni(i) = obj.params.ni_HTL;
                end
            end
            
            % 获取辐射复合系数
            if obj.position_dependent
                B = obj.getPositionDependentRadiativeCoefficient(x);
            else
                B = obj.B * ones(size(x));
            end
            
            % 计算辐射复合率
            R_rad = B .* (n .* p - ni.^2);
        end
        
        function R = addInterfaceRecombination(obj, R, n, p, x)
            % 添加界面复合
            
            % 使用已知的界面位置
            % 假设有两个界面：ETL/吸收层和吸收层/HTL
            interface_positions = [obj.params.L_ETL, obj.params.L_ETL + obj.params.L_absorber];
            
            % 对每个界面添加复合
            for i = 1:length(interface_positions)
                % 找到最接近界面的网格点
                [~, idx] = min(abs(x - interface_positions(i)));
                
                % 确定界面类型
                if i == 1  % ETL/吸收层界面
                    S_n = 1e3;  % 电子界面复合速率 [cm/s]
                    S_p = 5e2;  % 空穴界面复合速率 [cm/s]
                else  % 吸收层/HTL界面
                    S_n = 5e2;  % 电子界面复合速率 [cm/s]
                    S_p = 1e3;  % 空穴界面复合速率 [cm/s]
                end
                
                % 计算界面复合率
                R_interface = S_n * n(idx) + S_p * p(idx);
                
                % 将界面复合率添加到总复合率
                R(idx) = R(idx) + R_interface;
            end
        end
        
        function [tau_n, tau_p] = getPositionDependentSRHLifetimes(obj, x)
            % 获取位置依赖的SRH寿命
            
            % 初始化
            tau_n = obj.tau_n * ones(size(x));
            tau_p = obj.tau_p * ones(size(x));
            
            % 在不同区域设置不同的寿命
            idx_ETL = obj.params.idx_ETL;
            idx_abs = obj.params.idx_absorber;
            idx_HTL = obj.params.idx_HTL;
            
            % ETL中的寿命
            tau_n(idx_ETL) = 1e-9;  % 1 ns
            tau_p(idx_ETL) = 1e-9;
            
            % 吸收层中的寿命 (保持默认值)
            
            % HTL中的寿命
            tau_n(idx_HTL) = 1e-9;  % 1 ns
            tau_p(idx_HTL) = 1e-9;
            
            % 在界面附近降低寿命
            interfaces = obj.params.interfaces;
            for i = 1:length(interfaces)
                % 找到界面附近的点
                near_interface = abs(x - interfaces(i)) < 1e-6;  % 10 nm范围内
                
                % 降低界面附近的寿命
                tau_n(near_interface) = tau_n(near_interface) * 0.1;
                tau_p(near_interface) = tau_p(near_interface) * 0.1;
            end
        end
        
        function [Cn, Cp] = getPositionDependentAugerCoefficients(obj, x)
            % 获取位置依赖的Auger系数
            
            % 初始化
            Cn = obj.Cn * ones(size(x));
            Cp = obj.Cp * ones(size(x));
            
            % 在不同区域设置不同的系数
            idx_ETL = obj.params.idx_ETL;
            idx_abs = obj.params.idx_absorber;
            idx_HTL = obj.params.idx_HTL;
            
            % ETL中的系数
            Cn(idx_ETL) = 1e-29;
            Cp(idx_ETL) = 1e-29;
            
            % 吸收层中的系数 (保持默认值)
            
            % HTL中的系数
            Cn(idx_HTL) = 1e-29;
            Cp(idx_HTL) = 1e-29;
        end
        
        function B = getPositionDependentRadiativeCoefficient(obj, x)
            % 获取位置依赖的辐射复合系数
            
            % 初始化
            B = obj.B * ones(size(x));
            
            % 在不同区域设置不同的系数
            idx_ETL = obj.params.idx_ETL;
            idx_abs = obj.params.idx_absorber;
            idx_HTL = obj.params.idx_HTL;
            
            % ETL中的系数 (几乎没有辐射复合)
            B(idx_ETL) = 1e-12;
            
            % 吸收层中的系数 (保持默认值)
            
            % HTL中的系数 (几乎没有辐射复合)
            B(idx_HTL) = 1e-12;
        end
        
        function plotRecombinationRates(obj, n, p, x)
            % 绘制各种复合率
            
            % 计算各种复合率
            R_SRH = obj.calculateSRHRecombination(n, p, x);
            R_Auger = obj.calculateAugerRecombination(n, p, x);
            R_rad = obj.calculateRadiativeRecombination(n, p, x);
            R_total = R_SRH + R_Auger + R_rad;
            
            % 添加界面复合
            R_total = obj.addInterfaceRecombination(R_total, n, p, x);
            
            % 创建图形
            figure('Name', '复合率', 'Position', [100, 100, 800, 600]);
            
            % 绘制复合率
            semilogy(x*1e7, abs(R_SRH), 'r-', 'LineWidth', 2);
            hold on;
            semilogy(x*1e7, abs(R_Auger), 'g-', 'LineWidth', 2);
            semilogy(x*1e7, abs(R_rad), 'b-', 'LineWidth', 2);
            semilogy(x*1e7, abs(R_total), 'k-', 'LineWidth', 2);
            
            % 添加界面线
            for i = 1:length(obj.params.interfaces)
                x_int = obj.params.interfaces(i);
                plot([x_int, x_int]*1e7, [1e0, 1e25], 'k--');
            end
            
            % 添加层标签
            layers = {'ETL', '吸收层', 'HTL'};
            x_centers = zeros(1, 3);
            x_centers(1) = mean(x(obj.params.idx_ETL))*1e7;
            x_centers(2) = mean(x(obj.params.idx_absorber))*1e7;
            x_centers(3) = mean(x(obj.params.idx_HTL))*1e7;
            
            for i = 1:length(layers)
                text(x_centers(i), 1e20, layers{i}, ...
                     'HorizontalAlignment', 'center', 'FontWeight', 'bold');
            end
            
            % 设置图形属性
            xlabel('位置 (nm)');
            ylabel('复合率 (cm^{-3}s^{-1})');
            title('复合率分布');
            legend('SRH复合', 'Auger复合', '辐射复合', '总复合率', 'Location', 'best');
            grid on;
            
            % 设置坐标轴范围
            xlim([min(x), max(x)]*1e7);
            ylim([1e0, 1e25]);
        end
    end
end