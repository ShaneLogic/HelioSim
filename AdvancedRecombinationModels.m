classdef AdvancedRecombinationModels < handle
    % AdvancedRecombinationModels - 高级复合模型类
    % 包含SRH、Auger、辐射复合以及界面复合模型
    
    properties
        params      % 太阳能电池参数
        
        % SRH复合参数
        tau_n       % 电子寿命 [s]
        tau_p       % 空穴寿命 [s]
        Et          % 陷阱能级 [eV]
        
        % Auger复合参数
        Cn          % 电子Auger系数 [cm^6/s]
        Cp          % 空穴Auger系数 [cm^6/s]
        
        % 辐射复合参数
        B           % 辐射复合系数 [cm^3/s]
        
        % 界面复合参数
        S_n         % 电子界面复合速度 [cm/s]
        S_p         % 空穴界面复合速度 [cm/s]
    end
    
    methods
        function obj = AdvancedRecombinationModels(params)
            % 构造函数 - 初始化复合模型
            obj.params = params;
            
            % 设置默认SRH参数
            obj.tau_n = 1e-7;  % 电子寿命 [s]
            obj.tau_p = 1e-7;  % 空穴寿命 [s]
            obj.Et = 0;        % 中间带隙陷阱能级 [eV]
            
            % 设置默认Auger参数
            obj.Cn = 1e-30;    % 电子Auger系数 [cm^6/s]
            obj.Cp = 1e-30;    % 空穴Auger系数 [cm^6/s]
            
            % 设置默认辐射复合参数
            obj.B = 1e-10;     % 辐射复合系数 [cm^3/s]
            
            % 设置默认界面复合参数
            obj.S_n = 1e3;     % 电子界面复合速度 [cm/s]
            obj.S_p = 1e3;     % 空穴界面复合速度 [cm/s]
        end
        
        function setSRHParameters(obj, tau_n, tau_p, Et)
            % 设置SRH复合参数
            obj.tau_n = tau_n;
            obj.tau_p = tau_p;
            obj.Et = Et;
        end
        
        function setAugerParameters(obj, Cn, Cp)
            % 设置Auger复合参数
            obj.Cn = Cn;
            obj.Cp = Cp;
        end
        
        function setRadiativeParameters(obj, B)
            % 设置辐射复合参数
            obj.B = B;
        end
        
        function setInterfaceParameters(obj, S_n, S_p)
            % 设置界面复合参数
            obj.S_n = S_n;
            obj.S_p = S_p;
        end
        
        function R = calculateRecombinationRate(obj, n, p, x)
            % 计算总复合率
            % 输入：
            %   n - 电子密度 [cm^-3]
            %   p - 空穴密度 [cm^-3]
            %   x - 位置坐标 [cm]
            % 输出：
            %   R - 总复合率 [cm^-3 s^-1]
            
            % 计算位置相关的本征载流子浓度
            ni = zeros(size(x));
            for idx = 1:length(x)
                if x(idx) <= obj.params.L_ETL
                    % ETL区域
                    ni(idx) = sqrt(obj.params.Nc_ETL * obj.params.Nv_ETL) * ...
                        exp(-obj.params.Eg_ETL / (2 * obj.params.kb * obj.params.T / obj.params.q));
                elseif x(idx) <= obj.params.L_ETL + obj.params.L_absorber
                    % 吸收层区域
                    ni(idx) = sqrt(obj.params.Nc_abs * obj.params.Nv_abs) * ...
                        exp(-obj.params.Eg_abs / (2 * obj.params.kb * obj.params.T / obj.params.q));
                else
                    % HTL区域
                    ni(idx) = sqrt(obj.params.Nc_HTL * obj.params.Nv_HTL) * ...
                        exp(-obj.params.Eg_HTL / (2 * obj.params.kb * obj.params.T / obj.params.q));
                end
            end
            
            % 计算各种复合机制
            R_SRH = obj.calculateSRHRecombination(n, p, ni);
            R_Auger = obj.calculateAugerRecombination(n, p, ni);
            R_Radiative = obj.calculateRadiativeRecombination(n, p, ni);
            R_Interface = obj.calculateInterfaceRecombination(n, p, ni, x);
            
            % 计算总复合率
            R = R_SRH + R_Auger + R_Radiative + R_Interface;
        end
        
        function R_SRH = calculateSRHRecombination(obj, n, p, ni)
            % 计算SRH复合率
            % Shockley-Read-Hall复合模型
            
            % 计算陷阱能级
            Ei = zeros(size(ni));
            for i = 1:length(ni)
                Ei(i) = -obj.params.kb * obj.params.T / obj.params.q * log(ni(i));
            end
            
            % 计算SRH统计因子
            n1 = ni .* exp((obj.Et - Ei) / (obj.params.kb * obj.params.T / obj.params.q));
            p1 = ni .* exp(-(obj.Et - Ei) / (obj.params.kb * obj.params.T / obj.params.q));
            
            % 计算SRH复合率
            R_SRH = (n .* p - ni.^2) ./ (obj.tau_p .* (n + n1) + obj.tau_n .* (p + p1));
        end
        
        function R_Auger = calculateAugerRecombination(obj, n, p, ni)
            % 计算Auger复合率
            % Auger复合是三粒子过程
            
            R_Auger = obj.Cn .* n .* (n .* p - ni.^2) + obj.Cp .* p .* (n .* p - ni.^2);
        end
        
        function R_Radiative = calculateRadiativeRecombination(obj, n, p, ni)
            % 计算辐射复合率
            % 带间直接复合
            
            R_Radiative = obj.B .* (n .* p - ni.^2);
        end
        
        function R_Interface = calculateInterfaceRecombination(obj, n, p, ni, x)
            % 计算界面复合率
            % 在材料界面处的表面复合
            
            % 初始化界面复合率
            R_Interface = zeros(size(n));
            
            % 界面位置
            interface_positions = [obj.params.L_ETL, obj.params.L_ETL + obj.params.L_absorber];
            
            % 计算界面复合
            for interface_idx = 1:length(interface_positions)
                % 找到最接近界面的网格点
                [~, idx] = min(abs(x - interface_positions(interface_idx)));
                
                % 计算界面复合率
                % 使用表面复合速度模型
                R_Interface(idx) = (obj.S_n * obj.S_p) / (obj.S_n + obj.S_p) * (n(idx) * p(idx) - ni(idx)^2) / ni(idx);
            end
        end
    end
end
