% 物理常数
q = 1.6e-19;          % 电子电荷量 (C)
cmm = 1.0e6;          % 单位换算系数

% 读取电荷分布数据
data = readmatrix('rho.csv');
x = linspace(0, 1, size(data, 2)); % 定义 x 轴
time_steps = 1:10;                % 定义时间步长

% 定义颜色图，用于一致化曲线颜色
color_map = lines(length(time_steps));

% 初始化存储电荷分布函数、电势解和电场强度的变量
charge_density_fits = cell(1, 10);
potential_solutions = cell(1, 10);
electric_field_solutions = cell(1, 10);

% 使用 tiledlayout 创建三行图和一个独立图例
layout = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% 绘制电荷分布函数 f(x, t)
nexttile;
hold on;
for t = time_steps
    y = data(t, :) * q * cmm;       % 提取第 t 行数据
    f = chebfun.spline(x, y);       % 拟合电荷分布函数
    charge_density_fits{t} = f;     % 保存拟合结果
    plot(f, 'Color', color_map(t, :), 'LineWidth', 1); % 绘制曲线
end
title('Charge Density Distribution f(x, t)');
xlabel('x');
ylabel('f(x,t) (C/m^3)');
grid on;
set(gca, 'FontSize', 12);

% 求解 Poisson 方程并绘制电势解 ψ(x, t)
nexttile;
hold on;
for t = time_steps
    f = charge_density_fits{t}; % 直接读取拟合结果

    % 定义微分算子和作用域
    L = chebop(0, 1); % 定义作用域为 [0, 1]
    L.op = @(x, u) diff(u, 2); % 定义算子为二阶导数
    L.lbc = 0;  % 左边界条件：u(0) = 0
    L.rbc = 0;  % 右边界条件：u(1) = 0

    % 解方程 L[u] = f
    psi = L \ f; % 计算解 psi
    potential_solutions{t} = psi; % 保存电势解
    plot(psi, 'Color', color_map(t, :), 'LineWidth', 1); % 绘制曲线
end
title('Potential Functions \psi(x, t)');
xlabel('x');
ylabel('\psi(x, t) (V)');
grid on;
set(gca, 'FontSize', 12);

% 计算并绘制电场强度 E(x, t)
nexttile;
hold on;
for t = time_steps
    psi = potential_solutions{t}; % 读取保存的电势解
    E = -diff(psi);               % 计算电场强度（电势的一阶导）
    electric_field_solutions{t} = E; % 保存电场解
    plot(E, 'Color', color_map(t, :), 'LineWidth', 1); % 绘制曲线，无标记
end
title('Electric Field E(x, t)');
xlabel('x');
ylabel('E(x, t) (V/m)');
grid on;
set(gca, 'FontSize', 12);

% 添加统一图例，调整位置
lgd = legend(arrayfun(@(t) ['t = ' num2str(t)], time_steps, 'UniformOutput', false), ...
    'Orientation', 'vertical', 'FontSize', 12);
lgd.Layout.Tile = 'east'; % 放置图例到右侧
