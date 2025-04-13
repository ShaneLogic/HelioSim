classdef InterfaceHandler < handle
    % InterfaceHandler - 处理材料界面的载流子传输和复合
    
    properties
        params;          % 太阳能电池参数
        
        % 界面复合参数
        S_n_ETL_abs = 1e5;  % ETL/吸收层界面电子复合速率 [cm/s]
        S_p_ETL_abs = 1e5;  % ETL/吸收层界面空穴复合速率 [cm/s]
        S_n_abs_HTL = 1e5;  % 吸收层/HTL界面电子复合速率 [cm/s]
        S_p_abs_HTL = 1e5;  % 吸收层/HTL界面空穴复合速率 [cm/s]
    end
    
    methods
        function obj = InterfaceHandler(params)
            % 构造函数 - 初始化参数
            obj.params = params;
        end
        
        function [n, p] = applyInterfaceConditions(obj, n, p, phi)
            % 在界面处应用载流子连续性和电流连续性条件
            
            % 获取界面索引
            intf_idx = obj.params.idx_interfaces;
            
            % 处理ETL/吸收层界面
            if ~isempty(intf_idx) && length(intf_idx) >= 1
                idx_intf = intf_idx(1);
                
                % 获取界面两侧材料属性
                % 左侧(ETL)
                mu_n_L = obj.params.mu_n_ETL;
                mu_p_L = obj.params.mu_p_ETL;
                D_n_L = obj.params.D_n_ETL;
                D_p_L = obj.params.D_p_ETL;
                
                % 右侧(吸收层)
                mu_n_R = obj.params.mu_n_abs;
                mu_p_R = obj.params.mu_p_abs;
                D_n_R = obj.params.D_n_abs;
                D_p_R = obj.params.D_p_abs;
                
                % 使用热平衡界面模型处理电子
                % 计算界面两侧电子浓度差异系数
                [n(idx_intf-1), n(idx_intf+1)] = obj.calculateInterfaceConcentrations(n(idx_intf-1), n(idx_intf+1), ...
                    phi(idx_intf-1), phi(idx_intf+1), obj.params.chi_ETL, obj.params.chi_abs, 'electron');
                
                % 使用热平衡界面模型处理空穴
                [p(idx_intf-1), p(idx_intf+1)] = obj.calculateInterfaceConcentrations(p(idx_intf-1), p(idx_intf+1), ...
                    phi(idx_intf-1), phi(idx_intf+1), obj.params.chi_ETL + obj.params.Eg_ETL, obj.params.chi_abs + obj.params.Eg_abs, 'hole');
                
                % 界面处的浓度插值
                n(idx_intf) = sqrt(n(idx_intf-1) * n(idx_intf+1));
                p(idx_intf) = sqrt(p(idx_intf-1) * p(idx_intf+1));
            end
            
            % 处理吸收层/HTL界面
            if ~isempty(intf_idx) && length(intf_idx) >= 2
                idx_intf = intf_idx(2);
                
                % 获取界面两侧材料属性
                % 左侧(吸收层)
                mu_n_L = obj.params.mu_n_abs;
                mu_p_L = obj.params.mu_p_abs;
                D_n_L = obj.params.D_n_abs;
                D_p_L = obj.params.D_p_abs;
                
                % 右侧(HTL)
                mu_n_R = obj.params.mu_n_HTL;
                mu_p_R = obj.params.mu_p_HTL;
                D_n_R = obj.params.D_n_HTL;
                D_p_R = obj.params.D_p_HTL;
                
                % 使用热平衡界面模型处理电子
                [n(idx_intf-1), n(idx_intf+1)] = obj.calculateInterfaceConcentrations(n(idx_intf-1), n(idx_intf+1), ...
                    phi(idx_intf-1), phi(idx_intf+1), obj.params.chi_abs, obj.params.chi_HTL, 'electron');
                
                % 使用热平衡界面模型处理空穴
                [p(idx_intf-1), p(idx_intf+1)] = obj.calculateInterfaceConcentrations(p(idx_intf-1), p(idx_intf+1), ...
                    phi(idx_intf-1), phi(idx_intf+1), obj.params.chi_abs + obj.params.Eg_abs, obj.params.chi_HTL + obj.params.Eg_HTL, 'hole');
                
                % 界面处的浓度插值
                n(idx_intf) = sqrt(n(idx_intf-1) * n(idx_intf+1));
                p(idx_intf) = sqrt(p(idx_intf-1) * p(idx_intf+1));
            end
        end
        
        function [c_L, c_R] = calculateInterfaceConcentrations(obj, c_L_in, c_R_in, phi_L, phi_R, E_L, E_R, carrier_type)
            % 计算界面两侧的载流子浓度，考虑带阶的影响
            
            % 初始化输出
            c_L = c_L_in;
            c_R = c_R_in;
            
            % 计算带阶
            dE = E_R - E_L;
            
            % 将能级差转换为电压
            dV = dE / obj.params.q;
            
            % 计算电势差
            dPhi = phi_R - phi_L;
            
            % 计算有效势垒
            if strcmp(carrier_type, 'electron')
                barrier = dV - dPhi;
            else % 空穴
                barrier = -dV + dPhi;
            end
            
            % 应用界面条件
            % 如果载流子由左向右流动遇到下坡(barrier < 0)，则无阻碍
            % 如果载流子由左向右流动遇到上坡(barrier > 0)，则用玻尔兹曼因子调整
            if barrier > 0
                % 电子/空穴从左到右遇到上坡势垒
                c_R = c_L * exp(-barrier / obj.params.Vt);
            else
                % 电子/空穴从右到左遇到上坡势垒
                c_L = c_R * exp(barrier / obj.params.Vt);
            end
        end
        
        function R_interface = calculateInterfaceRecombination(obj, n, p, phi)
            % 计算界面复合率
            
            % 获取界面索引
            intf_idx = obj.params.idx_interfaces;
            Nx = length(n);
            
            % 初始化界面复合率
            R_interface = zeros(Nx, 1);
            
            % 处理ETL/吸收层界面
            if ~isempty(intf_idx) && length(intf_idx) >= 1
                idx_intf = intf_idx(1);
                
                % 获取界面两侧的载流子浓度
                n_L = n(idx_intf-1);
                n_R = n(idx_intf+1);
                p_L = p(idx_intf-1);
                p_R = p(idx_intf+1);
                
                % 界面处的浓度
                n_intf = n(idx_intf);
                p_intf = p(idx_intf);
                
                % 计算平衡态界面浓度
                ni_ETL = obj.params.ni_ETL;
                ni_abs = obj.params.ni_abs;
                ni_intf = sqrt(ni_ETL * ni_abs);
                
                % 计算界面复合率 (SRH模型)
                R_intf = (n_intf * p_intf - ni_intf^2) / ...
                         ((n_intf + ni_intf) / obj.S_p_ETL_abs + (p_intf + ni_intf) / obj.S_n_ETL_abs);
                
                % 将复合率赋值给界面点
                R_interface(idx_intf) = R_intf;
            end
            
            % 处理吸收层/HTL界面
            if ~isempty(intf_idx) && length(intf_idx) >= 2
                idx_intf = intf_idx(2);
                
                % 获取界面两侧的载流子浓度
                n_L = n(idx_intf-1);
                n_R = n(idx_intf+1);
                p_L = p(idx_intf-1);
                p_R = p(idx_intf+1);
                
                % 界面处的浓度
                n_intf = n(idx_intf);
                p_intf = p(idx_intf);
                
                % 计算平衡态界面浓度
                ni_abs = obj.params.ni_abs;
                ni_HTL = obj.params.ni_HTL;
                ni_intf = sqrt(ni_abs * ni_HTL);
                
                % 计算界面复合率 (SRH模型)
                R_intf = (n_intf * p_intf - ni_intf^2) / ...
                         ((n_intf + ni_intf) / obj.S_p_abs_HTL + (p_intf + ni_intf) / obj.S_n_abs_HTL);
                
                % 将复合率赋值给界面点
                R_interface(idx_intf) = R_intf;
            end
        end
    end
end