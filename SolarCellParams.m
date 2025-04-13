classdef SolarCellParams < handle
    % SolarCellParams - 太阳能电池参数类
    
    properties
        % 物理常数
        q = 1.602e-19;    % 基本电荷 [C]
        kb = 1.380e-23;   % 玻尔兹曼常数 [J/K]
        T = 300;          % 温度 [K]
        eps0 = 8.85e-14;  % 真空介电常数 [F/cm]
        Vt;               % 热电压 [V]
        
        % 设备几何参数
        L_ETL = 1e-5;     % ETL层厚度 [cm] (100 nm)
        L_absorber = 5e-5; % 吸收层厚度 [cm] (500 nm)
        L_HTL = 1e-5;     % HTL层厚度 [cm] (100 nm)
        L_total;          % 总厚度 [cm]
        
        % 网格参数
        Nx_ETL = 50;      % ETL层的网格点数
        Nx_absorber = 100; % 吸收层的网格点数
        Nx_HTL = 50;      % HTL层的网格点数
        Nx_total;         % 总网格点数
        x;                % 位置网格 [cm]
        dx;               % 网格间距 [cm]
        idx_ETL;          % ETL区域的索引
        idx_absorber;     % 吸收层区域的索引
        idx_HTL;          % HTL区域的索引
        idx_interfaces;   % 材料界面的索引
        
        % 材料参数 - ETL层 (如TiO2或ZnO)
        eps_ETL = 9.0;    % 相对介电常数
        chi_ETL = 4.0;    % 电子亲和能 [eV]
        Eg_ETL = 3.2;     % 带隙 [eV]
        Nc_ETL = 2.0e18;  % 导带有效态密度 [cm^-3]
        Nv_ETL = 1.8e19;  % 价带有效态密度 [cm^-3]
        mu_n_ETL = 100;   % 电子迁移率 [cm^2/Vs]
        mu_p_ETL = 25;    % 空穴迁移率 [cm^2/Vs]
        Nd_ETL = 1e17;    % 施主掺杂浓度 [cm^-3]
        Na_ETL = 0;       % 受主掺杂浓度 [cm^-3]
        ni_ETL;           % 本征载流子浓度 [cm^-3]
        tau_n_ETL = 1e-8; % 电子寿命 [s]
        tau_p_ETL = 1e-8; % 空穴寿命 [s]
        D_n_ETL;          % 电子扩散系数 [cm^2/s]
        D_p_ETL;          % 空穴扩散系数 [cm^2/s]
        Vbi_ETL;          % 内建电势 [V]
        
        % 材料参数 - 吸收层 (如钙钛矿或硅)
        eps_abs = 25;     % 相对介电常数
        chi_abs = 3.9;    % 电子亲和能 [eV]
        Eg_abs = 1.5;     % 带隙 [eV]
        Nc_abs = 1.0e18;  % 导带有效态密度 [cm^-3]
        Nv_abs = 1.0e18;  % 价带有效态密度 [cm^-3]
        mu_n_abs = 10;    % 电子迁移率 [cm^2/Vs]
        mu_p_abs = 10;    % 空穴迁移率 [cm^2/Vs]
        Nd_abs = 0;       % 施主掺杂浓度 [cm^-3]
        Na_abs = 0;       % 受主掺杂浓度 [cm^-3]
        ni_abs;           % 本征载流子浓度 [cm^-3]
        tau_n_abs = 1e-6; % 电子寿命 [s]
        tau_p_abs = 1e-6; % 空穴寿命 [s]
        D_n_abs;          % 电子扩散系数 [cm^2/s]
        D_p_abs;          % 空穴扩散系数 [cm^2/s]
        alpha_abs = 1e5;  % 吸收系数 [cm^-1]
        G_max = 1e21;     % 最大光生率 [cm^-3 s^-1]
        R_front = 0.05;   % 前表面反射率
        
        % 材料参数 - HTL层 (如Spiro-OMeTAD或PEDOT:PSS)
        eps_HTL = 3.0;    % 相对介电常数
        chi_HTL = 2.1;    % 电子亲和能 [eV]
        Eg_HTL = 3.0;     % 带隙 [eV]
        Nc_HTL = 2.5e18;  % 导带有效态密度 [cm^-3]
        Nv_HTL = 2.5e18;  % 价带有效态密度 [cm^-3]
        mu_n_HTL = 1;     % 电子迁移率 [cm^2/Vs]
        mu_p_HTL = 50;    % 空穴迁移率 [cm^2/Vs]
        Nd_HTL = 0;       % 施主掺杂浓度 [cm^-3]
        Na_HTL = 1e17;    % 受主掺杂浓度 [cm^-3]
        ni_HTL;           % 本征载流子浓度 [cm^-3]
        tau_n_HTL = 1e-8; % 电子寿命 [s]
        tau_p_HTL = 1e-8; % 空穴寿命 [s]
        D_n_HTL;          % 电子扩散系数 [cm^2/s]
        D_p_HTL;          % 空穴扩散系数 [cm^2/s]
        Vbi_HTL;          % 内建电势 [V]
        
        % 状态参数
        illumination = false; % 光照状态
        V_applied = 0;    % 施加电压 [V]
    end
    
    methods
        function obj = SolarCellParams()
            % 构造函数 - 初始化派生参数
            obj.Vt = obj.kb * obj.T / obj.q;
            obj.L_total = obj.L_ETL + obj.L_absorber + obj.L_HTL;
            
            % 设置计算网格
            obj.setupGrid();
            
            % 计算材料参数
            obj.setupMaterialParameters();
        end
        
        function setupGrid(obj)
            % 设置计算网格
            obj.Nx_total = obj.Nx_ETL + obj.Nx_absorber + obj.Nx_HTL - 2;
            
            % 为每层创建位置网格
            x_ETL = linspace(0, obj.L_ETL, obj.Nx_ETL);
            x_abs = linspace(obj.L_ETL, obj.L_ETL + obj.L_absorber, obj.Nx_absorber);
            x_HTL = linspace(obj.L_ETL + obj.L_absorber, obj.L_total, obj.Nx_HTL);
            
            % 合并网格，去除重复的界面点
            obj.x = [x_ETL(1:end-1), x_abs(1:end-1), x_HTL];
            obj.dx = diff(obj.x);
            
            % 存储每层的索引
            obj.idx_ETL = 1:(obj.Nx_ETL-1);
            obj.idx_absorber = obj.Nx_ETL:(obj.Nx_ETL+obj.Nx_absorber-2);
            obj.idx_HTL = (obj.Nx_ETL+obj.Nx_absorber-1):obj.Nx_total;
            
            % 存储界面索引
            obj.idx_interfaces = [obj.Nx_ETL, obj.Nx_ETL+obj.Nx_absorber-1];
        end
        
        function setupMaterialParameters(obj)
            % 计算派生材料参数
            
            % 本征载流子浓度
            obj.ni_ETL = sqrt(obj.Nc_ETL * obj.Nv_ETL) * exp(-obj.Eg_ETL/(2*obj.Vt));
            obj.ni_abs = sqrt(obj.Nc_abs * obj.Nv_abs) * exp(-obj.Eg_abs/(2*obj.Vt));
            obj.ni_HTL = sqrt(obj.Nc_HTL * obj.Nv_HTL) * exp(-obj.Eg_HTL/(2*obj.Vt));
            
            % 扩散系数 (爱因斯坦关系)
            obj.D_n_ETL = obj.mu_n_ETL * obj.Vt;
            obj.D_p_ETL = obj.mu_p_ETL * obj.Vt;
            obj.D_n_abs = obj.mu_n_abs * obj.Vt;
            obj.D_p_abs = obj.mu_p_abs * obj.Vt;
            obj.D_n_HTL = obj.mu_n_HTL * obj.Vt;
            obj.D_p_HTL = obj.mu_p_HTL * obj.Vt;
            
            % 计算内建电压
            obj.calculateBuiltInVoltage();
        end
        
        function calculateBuiltInVoltage(obj)
            % 计算内建电压
            % ETL侧内建电压 (ETL/吸收层界面)
            obj.Vbi_ETL = obj.Vt * log(obj.Nd_ETL * obj.Na_abs / (obj.ni_ETL * obj.ni_abs));
            if isnan(obj.Vbi_ETL) || obj.Vbi_ETL <= 0
                % 如果吸收层本征，使用不同的计算方式
                phi_f_ETL = obj.Vt * log(obj.Nd_ETL / obj.ni_ETL);
                phi_f_abs = 0; % 本征层费米能级在本征能级
                obj.Vbi_ETL = (obj.chi_ETL - obj.chi_abs)/obj.q + phi_f_ETL - phi_f_abs;
            end
            
            % HTL侧内建电压 (吸收层/HTL界面)
            obj.Vbi_HTL = obj.Vt * log(obj.Na_HTL * obj.Nd_abs / (obj.ni_HTL * obj.ni_abs));
            if isnan(obj.Vbi_HTL) || obj.Vbi_HTL <= 0
                % 如果吸收层本征，使用不同的计算方式
                phi_f_HTL = obj.Vt * log(obj.Na_HTL / obj.ni_HTL);
                phi_f_abs = 0; % 本征层费米能级在本征能级
                obj.Vbi_HTL = (obj.chi_abs - obj.chi_HTL)/obj.q + phi_f_abs - phi_f_HTL;
            end
        end
        
        function setAppliedVoltage(obj, voltage)
            % 设置施加电压
            obj.V_applied = voltage;
        end
        
        function setIllumination(obj, state)
            % 设置光照状态
            obj.illumination = state;
        end
        
        function setAbsorberThickness(obj, thickness)
            % 设置吸收层厚度 [cm]
            obj.L_absorber = thickness;
            obj.L_total = obj.L_ETL + obj.L_absorber + obj.L_HTL;
            obj.setupGrid();
        end
        
        function setAbsorberBandgap(obj, bandgap)
            % 设置吸收层带隙 [eV]
            obj.Eg_abs = bandgap;
            % 重新计算本征载流子浓度
            obj.ni_abs = sqrt(obj.Nc_abs * obj.Nv_abs) * exp(-obj.Eg_abs/(2*obj.Vt));
            % 更新内建电压
            obj.calculateBuiltInVoltage();
        end
        
        function setETLDoping(obj, doping)
            % 设置ETL层掺杂浓度 [cm^-3]
            obj.Nd_ETL = doping;
            % 更新内建电压
            obj.calculateBuiltInVoltage();
        end
        
        function setHTLDoping(obj, doping)
            % 设置HTL层掺杂浓度 [cm^-3]
            obj.Na_HTL = doping;
            % 更新内建电压
            obj.calculateBuiltInVoltage();
        end
        
        function setAbsorberMobility(obj, mobility)
            % 设置吸收层载流子迁移率 [cm^2/Vs]
            obj.mu_n_abs = mobility;
            obj.mu_p_abs = mobility;
            % 更新扩散系数
            obj.D_n_abs = obj.mu_n_abs * obj.Vt;
            obj.D_p_abs = obj.mu_p_abs * obj.Vt;
        end
    end
end