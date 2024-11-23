% 定义微分算子和边界条件
L = chebop(-1, 1);
L.op = @(x, u) diff(u, 2);
L.lbc = 0;  % 边界条件 u(-1) = 0
L.rbc = 1;  % 边界条件 u(1) = 0

% 定义参数
a_values = [0, 0.5, -0.5]; % a 的值
sigma2 = 0.01; % 高斯分布中的分母

% 循环计算不同 a 值的解
for i = 1:length(a_values)
    a = a_values(i); % 当前的 a 值
    f = chebfun(@(x) exp(-(x - a).^2 / sigma2)); % 定义 f(x)
    
    % 求解方程
    u = L \ f; % L[u] = f(x)
    
    % 绘制结果
    subplot(1, 3, i); % 创建子图
    plot(u, 'LineWidth', 1.5);
    title(['Solution for a = ', num2str(a)]);
    xlabel('x');
    ylabel('u(x)');
end
