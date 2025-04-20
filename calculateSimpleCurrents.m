function [Jn, Jp] = calculateSimpleCurrents(params, n, p, phi)
    % calculateSimpleCurrents - 计算电子和空穴的电流密度
    % 基于漫游-漫游方程的简化形式
    % 
    % 输入:
    %   params - 太阳能电池参数对象
    %   n - 电子密度 [cm^-3]
    %   p - 空穴密度 [cm^-3]
    %   phi - 电势 [V]
    %
    % 输出:
    %   Jn - 电子电流密度 [mA/cm^2]
    %   Jp - 空穴电流密度 [mA/cm^2]
    
    % 使用参数对象中的物理参数
    q = params.q; % 基本电荷 [C]
    kB = params.kb; % 玻尔兹曼常数 [J/K]
    T = params.T; % 温度 [K]
    
    % 初始化材料相关参数数组
    mu_n = zeros(size(params.x));
    mu_p = zeros(size(params.x));
    
    % 根据位置设置不同区域的迁移率
    for i = 1:length(params.x)
        if params.x(i) <= params.L_ETL
            % ETL区域
            mu_n(i) = params.mu_n_ETL;
            mu_p(i) = params.mu_p_ETL;
        elseif params.x(i) <= params.L_ETL + params.L_absorber
            % 吸收层区域
            mu_n(i) = params.mu_n_abs;
            mu_p(i) = params.mu_p_abs;
        else
            % HTL区域
            mu_n(i) = params.mu_n_HTL;
            mu_p(i) = params.mu_p_HTL;
        end
    end
    
    % 计算扬散系数 (爱因斯坦关系)
    Dn = kB * T / q * mu_n; % 电子扬散系数 [cm^2/s]
    Dp = kB * T / q * mu_p; % 空穴扬散系数 [cm^2/s]
    
    % 计算网格间距
    dx = mean(diff(params.x));
    
    % 初始化电流密度数组
    Jn = zeros(size(n));
    Jp = zeros(size(p));
    
    % 计算内部网格点的电流密度
    for idx = 2:length(n)-1
        % 使用中心差分计算电场
        E = -(phi(idx+1) - phi(idx-1)) / (2 * dx);
        
        % 计算电子电流密度 - 使用位置相关的迁移率
        drift_n = q * mu_n(idx) * n(idx) * E;
        diff_n = q * Dn(idx) * (n(idx+1) - n(idx-1)) / (2 * dx);
        Jn(idx) = drift_n - diff_n;
        
        % 计算空穴电流密度 - 使用位置相关的迁移率
        drift_p = q * mu_p(idx) * p(idx) * E;
        diff_p = q * Dp(idx) * (p(idx+1) - p(idx-1)) / (2 * dx);
        Jp(idx) = drift_p + diff_p;
    end
    
    % 边界处理
    Jn(1) = Jn(2);
    Jn(end) = Jn(end-1);
    Jp(1) = Jp(2);
    Jp(end) = Jp(end-1);
    
    % 对电流密度进行平滑处理，使用移动平均滤波
    window_size = 3;
    if length(Jn) > window_size
        Jn_smooth = zeros(size(Jn));
        Jp_smooth = zeros(size(Jp));
        
        for idx = 1:length(Jn)
            % 定义窗口范围
            start_idx = max(1, idx - floor(window_size/2));
            end_idx = min(length(Jn), idx + floor(window_size/2));
            
            % 计算窗口内的平均值
            Jn_smooth(idx) = mean(Jn(start_idx:end_idx));
            Jp_smooth(idx) = mean(Jp(start_idx:end_idx));
        end
        
        Jn = Jn_smooth;
        Jp = Jp_smooth;
    end
    
    % 限制电流密度的最大值以确保数值稳定性
    J_max = 20; % mA/cm^2
    Jn = max(min(Jn, J_max), -J_max);
    Jp = max(min(Jp, J_max), -J_max);
end
