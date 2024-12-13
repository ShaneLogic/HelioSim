% 定义变量
epsilon = 1.05e-12; % F/cm
q = 1.6e-19; % C
kT_over_q = 0.02585; % V
n0 = 1.45e10; % cm^-3
theta = 1e-16; % cm
m = 200; % grid points

% 初始化误差
tol = 1e-6; % 容差
error = Inf; % 初始误差设为无穷大

VP = -0.3475; % p型材料
VN = 0.3475;  % n型材料

% 定义掺杂浓度
NA = 10e-6;   % p型材料的受主掺杂浓度 (cm^-3)
ND = 10e-6;   % n型材料的施主掺杂浓度 (cm^-3)

% 初始化电势
V = zeros(m+2, 1);
V(1:floor(m/2)+1) = VP;
V(floor(m/2)+2:end) = VN;

% 预计算系数矩阵
Mj = 2 / theta^2 + (2 * q * n0 / (epsilon * kT_over_q)) * cosh(V(2:m+1) / kT_over_q);
CM = sparse(1:m, 1:m, Mj, m, m) + ...
     sparse(1:m-1, 2:m, (-1/theta^2) * ones(m-1, 1), m, m) + ...
     sparse(2:m, 1:m-1, (-1/theta^2) * ones(m-1, 1), m, m);

% 收敛性测试
errors = [];
while error > tol
    d2V_by_dx2 = (V(1:m) - 2*V(2:m+1) + V(3:m+2)) / theta^2; % 近似 (16)
    rho = q * (ND - NA - 2*n0 * sinh(V(2:m+1) / kT_over_q)); % 方程 (24)
    R = d2V_by_dx2 + rho / epsilon; % 方程 (21) 的右侧

    % 计算修正向量
    DV = CM \ R;
    % 更新电势向量
    V(2:m+1) = V(2:m+1) + DV;
    % 计算误差
    error = norm(DV, 2) / sqrt(m);
    errors = [errors, error];
end

% 可视化结果
figure;
subplot(2, 1, 1);
plot(V, '-o');
title('Potential Distribution');
xlabel('Grid Points');
ylabel('Potential (V)');

subplot(2, 1, 2);
semilogy(errors, '-x');
title('Convergence of Error');
xlabel('Iteration');
ylabel('Error');

subplot(2, 1, 1);