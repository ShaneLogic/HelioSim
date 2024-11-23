% 定义微分算子和作用域
L = chebop(-1, 1); % 定义作用域为 [-1, 1]
L.op = @(x, u) diff(u, 2); % 定义算子为二阶导数

% 设置边界条件
L.lbc = 0;  % 左边界条件：u(-1) = 0
L.rbc = 1;  % 右边界条件：u(1) = 1

% 定义方程右端函数为 Chebfun 对象
f = chebfun(@(x) sin(pi * x), [-1, 1]); % f(x) = sin(πx) (作为 Chebfun)

% 求解微分方程
u = L \ f; % 求解 L[u] = f

% 绘制解的图像
plot(u);
title('Solution of the Differential Equation'); % 标题
xlabel('x');           % x 轴标签
ylabel('u(x)');        % y 轴标签
