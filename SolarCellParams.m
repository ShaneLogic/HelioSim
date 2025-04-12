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
        
        % 材料参数、掺杂等其他属性...
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
        
        function setAppliedVoltage(obj, voltage)
            % 设置施加电压
            obj.V_applied = voltage;
        end
        
        function setIllumination(obj, state)
            % 设置光照状态
            obj.illumination = state;
        end
    end
end