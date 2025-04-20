classdef InterfaceHandlerOptimized < handle
    % InterfaceHandlerOptimized - 优化的界面处理类
    % 处理太阳能电池中的异质结界面，实现Scharfetter-Gummel离散化
    
    properties
        params          % 太阳能电池参数
        use_thermionic  % 是否使用热电子发射模型
    end
    
    methods
        function obj = InterfaceHandlerOptimized(params)
            % 构造函数 - 初始化界面处理器
            obj.params = params;
            
            % 默认使用热电子发射模型
            obj.use_thermionic = true;
            
            % 初始化界面位置
            % 如果参数中没有界面位置，则根据网格计算
            if ~isfield(params, 'idx_interfaces')
                params.idx_interfaces = [params.Nx_ETL, params.Nx_ETL + params.Nx_absorber];
            end
        end
        
        function [n_new, p_new] = applyInterfaceConditions(obj, n, p, phi)
            % 应用界面条件
            % n, p: 电子和空穴密度
            % phi: 电势
            % 返回应用界面条件后的载流子密度
            
            % 复制输入数组
            n_new = n;
            p_new = p;
            
            % 获取界面位置索引
            idx_interfaces = obj.params.idx_interfaces;
            
            % 对每个界面应用条件
            for i = 1:length(idx_interfaces)
                idx = idx_interfaces(i);
                
                % 确定界面类型
                if i == 1  % ETL/吸收层界面
                    % 获取界面两侧的材料参数
                    chi_left = obj.params.chi_ETL;
                    chi_right = obj.params.chi_abs;
                    Eg_left = obj.params.Eg_ETL;
                    Eg_right = obj.params.Eg_abs;
                    Nc_left = obj.params.Nc_ETL;
                    Nc_right = obj.params.Nc_abs;
                    Nv_left = obj.params.Nv_ETL;
                    Nv_right = obj.params.Nv_abs;
                    
                    % 界面复合速率
                    S_n = obj.params.S_n_ETL_abs;
                    S_p = obj.params.S_p_ETL_abs;
                else       % 吸收层/HTL界面
                    % 获取界面两侧的材料参数
                    chi_left = obj.params.chi_abs;
                    chi_right = obj.params.chi_HTL;
                    Eg_left = obj.params.Eg_abs;
                    Eg_right = obj.params.Eg_HTL;
                    Nc_left = obj.params.Nc_abs;
                    Nc_right = obj.params.Nc_HTL;
                    Nv_left = obj.params.Nv_abs;
                    Nv_right = obj.params.Nv_HTL;
                    
                    % 界面复合速率
                    S_n = obj.params.S_n_abs_HTL;
                    S_p = obj.params.S_p_abs_HTL;
                end
                
                % 计算能带不连续性
                dEc = obj.params.q_e * (chi_left - chi_right);  % 导带不连续性 [J]
                dEv = obj.params.q_e * ((chi_left + Eg_left) - (chi_right + Eg_right));  % 价带不连续性 [J]
                
                % 获取界面两侧的索引
                idx_left = idx - 1;
                idx_right = idx + 1;
                
                % 获取界面两侧的电势和载流子密度
                phi_left = phi(idx_left);
                phi_right = phi(idx_right);
                n_left = n(idx_left);
                n_right = n(idx_right);
                p_left = p(idx_left);
                p_right = p(idx_right);
                
                % 计算界面处的电场
                dx = obj.params.x(idx_right) - obj.params.x(idx_left);
                E_interface = -(phi_right - phi_left) / dx;  % [V/cm]
                
                if obj.use_thermionic
                    % 使用热电子发射模型
                    [n_left_new, n_right_new, p_left_new, p_right_new] = ...
                        obj.applyThermionicEmission(n_left, n_right, p_left, p_right, ...
                                                   phi_left, phi_right, dEc, dEv, ...
                                                   Nc_left, Nc_right, Nv_left, Nv_right);
                else
                    % 使用简单的连续性条件
                    [n_left_new, n_right_new, p_left_new, p_right_new] = ...
                        obj.applyContinuityConditions(n_left, n_right, p_left, p_right, ...
                                                     phi_left, phi_right, dEc, dEv, ...
                                                     Nc_left, Nc_right, Nv_left, Nv_right);
                end
                
                % 更新界面两侧的载流子密度
                n_new(idx_left) = n_left_new;
                n_new(idx_right) = n_right_new;
                p_new(idx_left) = p_left_new;
                p_new(idx_right) = p_right_new;
                
                % 在界面处插值
                n_new(idx) = sqrt(n_left_new * n_right_new);
                p_new(idx) = sqrt(p_left_new * p_right_new);
            end
        end
        
        function [n_left, n_right, p_left, p_right] = applyThermionicEmission(obj, ...
                n_left_old, n_right_old, p_left_old, p_right_old, ...
                phi_left, phi_right, dEc, dEv, ...
                Nc_left, Nc_right, Nv_left, Nv_right)
            % 应用热电子发射模型
            % 考虑能带不连续性和热电子发射电流
            
            % 热电压
            Vt = obj.params.kb * obj.params.T / obj.params.q_e;
            
            % 计算界面处的电势差
            dphi = phi_right - phi_left;
            
            % 计算热电子发射电流
            % 电子电流
            Jn_left_to_right = obj.params.q_e * obj.params.vth_n * n_left_old * exp(-dEc/(obj.params.kb*obj.params.T));
            Jn_right_to_left = obj.params.q_e * obj.params.vth_n * n_right_old;
            
            % 空穴电流
            Jp_left_to_right = obj.params.q_e * obj.params.vth_p * p_left_old;
            Jp_right_to_left = obj.params.q_e * obj.params.vth_p * p_right_old * exp(-dEv/(obj.params.kb*obj.params.T));
            
            % 净电流
            Jn_net = Jn_left_to_right - Jn_right_to_left;
            Jp_net = Jp_left_to_right - Jp_right_to_left;
            
            % 更新载流子密度
            % 使用指数形式的Scharfetter-Gummel离散化
            B_plus = @(x) x ./ (exp(x) - 1);
            B_minus = @(x) B_plus(-x);
            
            % 电子密度
            n_left = n_left_old;
            n_right = n_right_old;
            
            % 空穴密度
            p_left = p_left_old;
            p_right = p_right_old;
            
            % 考虑能带不连续性
            % 电子准费米能级连续
            Efn_left = obj.params.kb * obj.params.T * log(n_left_old / Nc_left);
            Efn_right = Efn_left;
            n_right = Nc_right * exp(Efn_right / (obj.params.kb * obj.params.T));
            
            % 空穴准费米能级连续
            Efp_left = -obj.params.kb * obj.params.T * log(p_left_old / Nv_left);
            Efp_right = Efp_left;
            p_right = Nv_right * exp(-Efp_right / (obj.params.kb * obj.params.T));
        end
        
        function [n_left, n_right, p_left, p_right] = applyContinuityConditions(obj, ...
                n_left_old, n_right_old, p_left_old, p_right_old, ...
                phi_left, phi_right, dEc, dEv, ...
                Nc_left, Nc_right, Nv_left, Nv_right)
            % 应用简单的连续性条件
            % 假设准费米能级在界面处连续
            
            % 热电压
            Vt = obj.params.kb * obj.params.T / obj.params.q_e;
            
            % 电子准费米能级连续
            Efn_left = obj.params.kb * obj.params.T * log(n_left_old / Nc_left);
            Efn_right = Efn_left;
            n_right = Nc_right * exp(Efn_right / (obj.params.kb * obj.params.T));
            n_left = n_left_old;
            
            % 空穴准费米能级连续
            Efp_left = -obj.params.kb * obj.params.T * log(p_left_old / Nv_left);
            Efp_right = Efp_left;
            p_right = Nv_right * exp(-Efp_right / (obj.params.kb * obj.params.T));
            p_left = p_left_old;
        end
        
        function [Jn, Jp] = calculateInterfaceCurrents(obj, n, p, phi)
            % 计算界面处的电流密度
            % 使用Scharfetter-Gummel离散化
            
            % 获取界面位置索引
            idx_interfaces = obj.params.idx_interfaces;
            
            % 初始化电流数组
            Jn = zeros(size(n));
            Jp = zeros(size(p));
            
            % 对每个界面计算电流
            for i = 1:length(idx_interfaces)
                idx = idx_interfaces(i);
                
                % 确定界面类型
                if i == 1  % ETL/吸收层界面
                    % 获取界面两侧的材料参数
                    mu_n_left = obj.params.mu_n_ETL;
                    mu_n_right = obj.params.mu_n_abs;
                    mu_p_left = obj.params.mu_p_ETL;
                    mu_p_right = obj.params.mu_p_abs;
                    
                    % 能带不连续性
                    dEc = obj.params.q_e * (obj.params.chi_ETL - obj.params.chi_abs);
                    dEv = obj.params.q_e * ((obj.params.chi_ETL + obj.params.Eg_ETL) - ...
                                           (obj.params.chi_abs + obj.params.Eg_abs));
                else       % 吸收层/HTL界面
                    % 获取界面两侧的材料参数
                    mu_n_left = obj.params.mu_n_abs;
                    mu_n_right = obj.params.mu_n_HTL;
                    mu_p_left = obj.params.mu_p_abs;
                    mu_p_right = obj.params.mu_p_HTL;
                    
                    % 能带不连续性
                    dEc = obj.params.q_e * (obj.params.chi_abs - obj.params.chi_HTL);
                    dEv = obj.params.q_e * ((obj.params.chi_abs + obj.params.Eg_abs) - ...
                                           (obj.params.chi_HTL + obj.params.Eg_HTL));
                end
                
                % 获取界面两侧的索引
                idx_left = idx - 1;
                idx_right = idx + 1;
                
                % 获取界面两侧的电势和载流子密度
                phi_left = phi(idx_left);
                phi_right = phi(idx_right);
                n_left = n(idx_left);
                n_right = n(idx_right);
                p_left = p(idx_left);
                p_right = p(idx_right);
                
                % 计算界面处的电场
                dx = obj.params.x(idx_right) - obj.params.x(idx_left);
                E_interface = -(phi_right - phi_left) / dx;  % [V/cm]
                
                % 使用Scharfetter-Gummel离散化计算电流
                % 热电压
                Vt = obj.params.kb * obj.params.T / obj.params.q_e;
                
                % 电势差（包括能带不连续性）
                dphi_n = (phi_right - phi_left) + dEc/obj.params.q_e;
                dphi_p = (phi_right - phi_left) - dEv/obj.params.q_e;
                
                % Bernoulli函数
                B_plus = @(x) x ./ (exp(x) - 1);
                B_minus = @(x) B_plus(-x);
                
                % 平均迁移率
                mu_n_avg = 2 * mu_n_left * mu_n_right / (mu_n_left + mu_n_right);
                mu_p_avg = 2 * mu_p_left * mu_p_right / (mu_p_left + mu_p_right);
                
                % 计算电流密度
                Jn(idx) = obj.params.q_e * mu_n_avg * Vt / dx * ...
                          (n_left * B_minus(dphi_n/Vt) - n_right * B_plus(dphi_n/Vt));
                
                Jp(idx) = obj.params.q_e * mu_p_avg * Vt / dx * ...
                          (p_right * B_minus(dphi_p/Vt) - p_left * B_plus(dphi_p/Vt));
            end
        end
        
        function enableThermionicEmission(obj, enable)
            % 启用或禁用热电子发射模型
            obj.use_thermionic = enable;
        end
    end
end