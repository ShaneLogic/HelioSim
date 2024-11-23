% 物理常数
q = 1.6e-19; % 电子电荷量
epsilon = 8.85e-12; % 介电常数
mu_n = 0.1; % 电子迁移率
mu_p = 0.05; % 空穴迁移率
D_n = 0.025; % 电子扩散系数
D_p = 0.015; % 空穴扩散系数

% 定义Chebfun操作区间
a = 0;
b = 1;
L = chebop(a,b);

% 定义变量
x = chebfun('x',[a,b]);

% 初始猜测函数（可根据问题调整）
varphi_0 = chebfun(@(x) 0*x, [a,b]);
n_0 = chebfun(@(x) 1e10*ones(size(x)), [a,b]);
p_0 = chebfun(@(x) 1e10*ones(size(x)), [a,b]);
E_0 = chebfun(@(x) 0*x, [a,b]); % 新的未知函数 E

% 定义耦合方程
rho = q*(p_0 - n_0); % 电荷密度

L.op = @(x,varphi,E,n,p) [diff(varphi) + E;...
    diff(E) + rho/epsilon;...
    diff(q*mu_n*n.*E + q*D_n*diff(n));...
    diff(q*mu_p*p.*E - q*D_p*diff(p))];

% 设置边界条件（示例：Dirichlet边界条件）
varphi_a = 0;
varphi_b = 1;
n_a = 1e10;
n_b = 1e10;
p_a = 1e10;
p_b = 1e10;

L.lbc = @(varphi,E,n,p) [varphi - varphi_a; n - n_a; p - p_a];
L.rbc = @(varphi,E,n,p) [varphi - varphi_b; n - n_b; p - p_b];

% 求解耦合方程
[varphi,E,n,p] = L\[0;0;0;0];

% 绘图
subplot(4,1,1)
plot(varphi)
title('Electric Potential \varphi(x)')
xlabel('x')
ylabel('V')

subplot(4,1,2)
plot(E)
title('Electric Field E(x)')
xlabel('x')
ylabel('V/m')

subplot(4,1,3)
plot(n)
title('Electron Concentration n(x)')
xlabel('x')
ylabel('m^{-3}')

subplot(4,1,4)
plot(p)
title('Hole Concentration p(x)')
xlabel('x')
ylabel('m^{-3}')