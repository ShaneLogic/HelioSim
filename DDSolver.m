% DDSolver.m
classdef DDSolver < handle
    properties
        params;         % 太阳能电池参数
        config;         % 模拟配置
        optical_gen;    % 光生载流子生成对象
        recomb_models;  % 复合模型对象
        interface;      % 界面处理对象
    end
    
    methods
        function obj = DDSolver(params, config)
            % 构造函数 - 初始化求解器
            obj.params = params;
            obj.config = config;
            
            % 初始化其他子模块
            obj.optical_gen = OpticalGeneration(params);
            obj.recomb_models = RecombinationModels(params);
            obj.interface = InterfaceHandler(params);
        end
        
        function results = solve(obj)
            % 求解漂移-扩散方程组
            % 首先计算初始条件，然后迭代求解至稳态
            
            % 获取网格和参数
            x = obj.params.x;
            Nx = length(x);
            
            % 时间步长设置
            t_start = obj.config.t_start;
            t_end = obj.config.t_end;
            Nt = obj.config.num_time_steps;
            t = linspace(t_start, t_end, Nt);
            dt = diff(t);
            
            % 初始化结果数组
            results = struct();
            results.t = t;
            results.phi = zeros(Nt, Nx);
            results.n = zeros(Nt, Nx);
            results.p = zeros(Nt, Nx);
            results.Jn = zeros(Nt, Nx);
            results.Jp = zeros(Nt, Nx);
            results.J_total = zeros(Nt, Nx);
            
            % 计算平衡态初始条件
            [phi0, n0, p0] = obj.calculateEquilibrium();
            results.phi(1,:) = phi0;
            results.n(1,:) = n0;
            results.p(1,:) = p0;
            
            % 计算光生载流子生成率
            G = obj.optical_gen.calculateGeneration();
            
            % 迭代求解时间发展
            for i = 2:Nt
                % 打印进度
                if mod(i, 10) == 0
                    fprintf('时间步 %d/%d 完成, t = %.2e s\n', i, Nt, t(i));
                end
                
                % 上一个时间步的解
                phi_prev = results.phi(i-1,:)';
                n_prev = results.n(i-1,:)';
                p_prev = results.p(i-1,:)';
                
                % 计算复合率
                R = obj.recomb_models.calculateRecombination(n_prev, p_prev);
                
                % 求解泊松方程
                phi_new = obj.solvePoissonEquation(phi_prev, n_prev, p_prev);
                
                % 求解电子连续性方程
                n_new = obj.solveElectronContinuity(n_prev, p_prev, phi_new, G, R, dt(i-1));
                
                % 求解空穴连续性方程
                p_new = obj.solveHoleContinuity(n_prev, p_prev, phi_new, G, R, dt(i-1));
                
                % 存储结果
                results.phi(i,:) = phi_new';
                results.n(i,:) = n_new';
                results.p(i,:) = p_new';
                
                % 计算电流密度
                [Jn, Jp] = obj.calculateCurrentDensities(n_new, p_new, phi_new);
                results.Jn(i,:) = Jn';
                results.Jp(i,:) = Jp';
                results.J_total(i,:) = (Jn + Jp)';
                
                % 检查收敛性
                if i > 10 && obj.isConverged(results, i, obj.config.rel_tol, obj.config.abs_tol)
                    fprintf('求解在 t = %.2e s 收敛\n', t(i));
                    % 将最后一个解复制到剩余时间步
                    for j = i+1:Nt
                        results.phi(j,:) = results.phi(i,:);
                        results.n(j,:) = results.n(i,:);
                        results.p(j,:) = results.p(i,:);
                        results.Jn(j,:) = results.Jn(i,:);
                        results.Jp(j,:) = results.Jp(i,:);
                        results.J_total(j,:) = results.J_total(i,:);
                    end
                    break;
                end
            end
            
            % 计算额外的结果数据
            results = obj.calculateAdditionalResults(results);
        end
        
        function [phi0, n0, p0] = calculateEquilibrium(obj)
            % 计算平衡态的电势和载流子浓度
            
            % 获取网格和参数
            x = obj.params.x;
            Nx = length(x);
            
            % 初始化电势和载流子浓度
            phi0 = zeros(Nx, 1);
            n0 = zeros(Nx, 1);
            p0 = zeros(Nx, 1);
            
            % 在各区域填充平衡载流子浓度
            % ETL区域
            idx = obj.params.idx_ETL;
            n0(idx) = obj.params.Nd_ETL;
            p0(idx) = obj.params.ni_ETL^2 ./ n0(idx);
            
            % 吸收层区域
            idx = obj.params.idx_absorber;
            n0(idx) = obj.params.Nd_abs;
            p0(idx) = obj.params.Na_abs;
            if obj.params.Nd_abs == 0 && obj.params.Na_abs == 0
                n0(idx) = obj.params.ni_abs;
                p0(idx) = obj.params.ni_abs;
            elseif obj.params.Nd_abs > 0
                n0(idx) = obj.params.Nd_abs;
                p0(idx) = obj.params.ni_abs^2 ./ n0(idx);
            else
                p0(idx) = obj.params.Na_abs;
                n0(idx) = obj.params.ni_abs^2 ./ p0(idx);
            end
            
            % HTL区域
            idx = obj.params.idx_HTL;
            p0(idx) = obj.params.Na_HTL;
            n0(idx) = obj.params.ni_HTL^2 ./ p0(idx);
            
            % 求解平衡态的泊松方程
            phi0 = obj.solvePoissonEquation(zeros(Nx, 1), n0, p0);
            
            % 处理界面
            [n0, p0] = obj.interface.applyInterfaceConditions(n0, p0, phi0);
        end
        
        function phi = solvePoissonEquation(obj, phi_init, n, p)
            % 数值求解泊松方程: ∇²φ = q/ε₀ε * (p - n + N⁺ - N⁻)
            
            % 获取网格和材料参数
            x = obj.params.x;
            Nx = length(x);
            dx = diff(x);
            
            % 获取掺杂浓度
            Nd = zeros(Nx, 1);
            Na = zeros(Nx, 1);
            
            % ETL区域
            idx = obj.params.idx_ETL;
            Nd(idx) = obj.params.Nd_ETL;
            
            % 吸收层区域
            idx = obj.params.idx_absorber;
            Nd(idx) = obj.params.Nd_abs;
            Na(idx) = obj.params.Na_abs;
            
            % HTL区域
            idx = obj.params.idx_HTL;
            Na(idx) = obj.params.Na_HTL;
            
            % 介电常数
            eps = zeros(Nx, 1);
            eps(obj.params.idx_ETL) = obj.params.eps_ETL;
            eps(obj.params.idx_absorber) = obj.params.eps_abs;
            eps(obj.params.idx_HTL) = obj.params.eps_HTL;
            
            % 构建系数矩阵
            A = spdiags(ones(Nx,1) * [-1 2 -1], [-1 0 1], Nx, Nx);
            
            % 处理边界条件
            % 左边界: 固定电势
            A(1,1) = 1;
            A(1,2) = 0;
            
            % 右边界: 固定电势
            A(Nx,Nx) = 1;
            A(Nx,Nx-1) = 0;
            
            % 调整非均匀网格的系数
            for i = 2:Nx-1
                dx_left = x(i) - x(i-1);
                dx_right = x(i+1) - x(i);
                dx_avg = (dx_left + dx_right) / 2;
                
                A(i,i-1) = -2 / (dx_left * (dx_left + dx_right));
                A(i,i) = 2 / (dx_left * dx_right);
                A(i,i+1) = -2 / (dx_right * (dx_left + dx_right));
            end
            
            % 构建右侧向量
            b = zeros(Nx, 1);
            
            % 内部节点:泊松方程右侧
            for i = 2:Nx-1
                b(i) = obj.params.q / obj.params.eps0 * ...
                        (p(i) - n(i) + Nd(i) - Na(i)) / eps(i);
            end
            
            % 边界条件
            phi_left = obj.params.Vbi_ETL;  % 左边界电势
            phi_right = -obj.params.Vbi_HTL + obj.params.V_applied;  % 右边界电势
            
            b(1) = phi_left;
            b(Nx) = phi_right;
            
            % 求解线性方程组
            phi = A \ b;
        end
        
        function n_new = solveElectronContinuity(obj, n, p, phi, G, R, dt)
            % 求解电子连续性方程: ∂n/∂t = ∇·(Dn∇n + μn·n·∇φ) + G - R
            
            % 获取网格和参数
            x = obj.params.x;
            Nx = length(x);
            
            % 提取电子迁移率和扩散系数
            mu_n = zeros(Nx, 1);
            mu_n(obj.params.idx_ETL) = obj.params.mu_n_ETL;
            mu_n(obj.params.idx_absorber) = obj.params.mu_n_abs;
            mu_n(obj.params.idx_HTL) = obj.params.mu_n_HTL;
            
            D_n = zeros(Nx, 1);
            D_n(obj.params.idx_ETL) = obj.params.D_n_ETL;
            D_n(obj.params.idx_absorber) = obj.params.D_n_abs;
            D_n(obj.params.idx_HTL) = obj.params.D_n_HTL;
            
            % 构建隐式方法系数矩阵
            A = spdiags(ones(Nx,1) * [-1 2 -1], [-1 0 1], Nx, Nx);
            
            % 调整内部节点的系数
            for i = 2:Nx-1
                dx_left = x(i) - x(i-1);
                dx_right = x(i+1) - x(i);
                dx_avg = (dx_left + dx_right) / 2;
                
                % 扩散项系数
                A(i,i-1) = -D_n(i) / (dx_left * dx_avg);
                A(i,i) = D_n(i) * (1/dx_left + 1/dx_right) / dx_avg;
                A(i,i+1) = -D_n(i) / (dx_right * dx_avg);
                
                % 漂移项
                dphi_dx_left = (phi(i) - phi(i-1)) / dx_left;
                dphi_dx_right = (phi(i+1) - phi(i)) / dx_right;
                
                A(i,i-1) = A(i,i-1) + mu_n(i) * dphi_dx_left / dx_avg;
                A(i,i+1) = A(i,i+1) - mu_n(i) * dphi_dx_right / dx_avg;
            end
            
            % 边界条件
            % 左边界: 表面复合
            A(1,1) = 1;
            A(1,2) = 0;
            
            % 右边界: 表面复合
            A(Nx,Nx) = 1;
            A(Nx,Nx-1) = 0;
            
            % 构建右侧向量
            b = n / dt + G - R;
            
            % 边界条件
            b(1) = n(1);  % 左边界
            b(Nx) = 0;    % 右边界(极低电子浓度)
            
            % 添加网格系数
            A = spdiags(ones(Nx,1) / dt, 0, Nx, Nx) + A;
            
            % 求解线性方程组
            n_new = A \ b;
            
            % 处理界面
            [n_new, ~] = obj.interface.applyInterfaceConditions(n_new, p, phi);
        end
        
        function p_new = solveHoleContinuity(obj, n, p, phi, G, R, dt)
            % 求解空穴连续性方程: ∂p/∂t = ∇·(Dp∇p - μp·p·∇φ) + G - R
            
            % 获取网格和参数
            x = obj.params.x;
            Nx = length(x);
            
            % 提取空穴迁移率和扩散系数
            mu_p = zeros(Nx, 1);
            mu_p(obj.params.idx_ETL) = obj.params.mu_p_ETL;
            mu_p(obj.params.idx_absorber) = obj.params.mu_p_abs;
            mu_p(obj.params.idx_HTL) = obj.params.mu_p_HTL;
            
            D_p = zeros(Nx, 1);
            D_p(obj.params.idx_ETL) = obj.params.D_p_ETL;
            D_p(obj.params.idx_absorber) = obj.params.D_p_abs;
            D_p(obj.params.idx_HTL) = obj.params.D_p_HTL;
            
            % 构建隐式方法系数矩阵
            A = spdiags(ones(Nx,1) * [-1 2 -1], [-1 0 1], Nx, Nx);
            
            % 调整内部节点的系数
            for i = 2:Nx-1
                dx_left = x(i) - x(i-1);
                dx_right = x(i+1) - x(i);
                dx_avg = (dx_left + dx_right) / 2;
                
                % 扩散项系数
                A(i,i-1) = -D_p(i) / (dx_left * dx_avg);
                A(i,i) = D_p(i) * (1/dx_left + 1/dx_right) / dx_avg;
                A(i,i+1) = -D_p(i) / (dx_right * dx_avg);
                
                % 漂移项
                dphi_dx_left = (phi(i) - phi(i-1)) / dx_left;
                dphi_dx_right = (phi(i+1) - phi(i)) / dx_right;
                
                A(i,i-1) = A(i,i-1) - mu_p(i) * dphi_dx_left / dx_avg;
                A(i,i+1) = A(i,i+1) + mu_p(i) * dphi_dx_right / dx_avg;
            end
            
            % 边界条件
            % 左边界: 表面复合
            A(1,1) = 1;
            A(1,2) = 0;
            
            % 右边界: 表面复合
            A(Nx,Nx) = 1;
            A(Nx,Nx-1) = 0;
            
            % 构建右侧向量
            b = p / dt + G - R;
            
            % 边界条件
            b(1) = 0;     % 左边界(极低空穴浓度)
            b(Nx) = p(Nx); % 右边界
            
            % 添加网格系数
            A = spdiags(ones(Nx,1) / dt, 0, Nx, Nx) + A;
            
            % 求解线性方程组
            p_new = A \ b;
            
            % 处理界面
            [~, p_new] = obj.interface.applyInterfaceConditions(n, p_new, phi);
        end
        
        function [Jn, Jp] = calculateCurrentDensities(obj, n, p, phi)
            % 计算电子和空穴的电流密度
            
            % 获取网格和参数
            x = obj.params.x;
            Nx = length(x);
            
            % 初始化电流密度数组
            Jn = zeros(Nx, 1);
            Jp = zeros(Nx, 1);
            
            % 提取迁移率
            mu_n = zeros(Nx, 1);
            mu_n(obj.params.idx_ETL) = obj.params.mu_n_ETL;
            mu_n(obj.params.idx_absorber) = obj.params.mu_n_abs;
            mu_n(obj.params.idx_HTL) = obj.params.mu_n_HTL;
            
            mu_p = zeros(Nx, 1);
            mu_p(obj.params.idx_ETL) = obj.params.mu_p_ETL;
            mu_p(obj.params.idx_absorber) = obj.params.mu_p_abs;
            mu_p(obj.params.idx_HTL) = obj.params.mu_p_HTL;
            
            % 提取扩散系数
            D_n = zeros(Nx, 1);
            D_n(obj.params.idx_ETL) = obj.params.D_n_ETL;
            D_n(obj.params.idx_absorber) = obj.params.D_n_abs;
            D_n(obj.params.idx_HTL) = obj.params.D_n_HTL;
            
            D_p = zeros(Nx, 1);
            D_p(obj.params.idx_ETL) = obj.params.D_p_ETL;
            D_p(obj.params.idx_absorber) = obj.params.D_p_abs;
            D_p(obj.params.idx_HTL) = obj.params.D_p_HTL;
            
            % 计算内部节点的电流密度
            for i = 2:Nx-1
                dx_left = x(i) - x(i-1);
                dx_right = x(i+1) - x(i);
                
                % 电势梯度
                dphi_dx = (phi(i+1) - phi(i-1)) / (dx_left + dx_right);
                
                % 载流子密度梯度
                dn_dx = (n(i+1) - n(i-1)) / (dx_left + dx_right);
                dp_dx = (p(i+1) - p(i-1)) / (dx_left + dx_right);
                
                % 电子电流密度: J_n = q·μ_n·n·E + q·D_n·∇n
                Jn(i) = obj.params.q * (mu_n(i) * n(i) * (-dphi_dx) + D_n(i) * dn_dx);
                
                % 空穴电流密度: J_p = q·μ_p·p·E - q·D_p·∇p
                Jp(i) = obj.params.q * (mu_p(i) * p(i) * dphi_dx - D_p(i) * dp_dx);
            end
            
            % 处理边界
            Jn(1) = Jn(2);
            Jn(Nx) = Jn(Nx-1);
            
            Jp(1) = Jp(2);
            Jp(Nx) = Jp(Nx-1);
        end
        
        function results = calculateAdditionalResults(obj, results)
            % 计算附加结果: 能带图、电场等
            
            % 获取网格和参数
            x = obj.params.x;
            Nx = length(x);
            
            % 获取最后时间步的结果
            phi = results.phi(end,:)';
            n = results.n(end,:)';
            p = results.p(end,:)';
            
            % 计算接合区能带弯曲
            % 电子亲和能
            chi = zeros(Nx, 1);
            chi(obj.params.idx_ETL) = obj.params.chi_ETL;
            chi(obj.params.idx_absorber) = obj.params.chi_abs;
            chi(obj.params.idx_HTL) = obj.params.chi_HTL;
            
            % 能带计算
            Ec = -obj.params.q * phi + chi;      % 导带底
            Ev = Ec - zeros(Nx, 1);              % 价带顶
            
            % 填充各区域的带隙
            Ev(obj.params.idx_ETL) = Ec(obj.params.idx_ETL) - obj.params.Eg_ETL;
            Ev(obj.params.idx_absorber) = Ec(obj.params.idx_absorber) - obj.params.Eg_abs;
            Ev(obj.params.idx_HTL) = Ec(obj.params.idx_HTL) - obj.params.Eg_HTL;
            
            % 计算费米能级
            Ef = obj.calculateFermiLevel(n, p);
            
            % 存储能带结果
            results.Ec = Ec;
            results.Ev = Ev;
            results.Ef = Ef;
            
            % 计算电场
            E_field = zeros(Nx, 1);
            for i = 2:Nx-1
                dx_left = x(i) - x(i-1);
                dx_right = x(i+1) - x(i);
                E_field(i) = -(phi(i+1) - phi(i-1)) / (dx_left + dx_right);
            end
            E_field(1) = E_field(2);
            E_field(Nx) = E_field(Nx-1);
            
            results.E_field = E_field;
        end
        
        function Ef = calculateFermiLevel(obj, n, p)
            % 计算费米能级
            
            % 获取网格和参数
            Nx = length(obj.params.x);
            
            % 初始化费米能级
            Ef = zeros(Nx, 1);
            
            % 计算各区域的费米能级
            % ETL区域
            idx = obj.params.idx_ETL;
            Ec_ETL = obj.params.chi_ETL;
            Nc_ETL = obj.params.Nc_ETL;
            Ef(idx) = Ec_ETL + obj.params.Vt * log(n(idx) ./ Nc_ETL);
            
            % 吸收层区域
            idx = obj.params.idx_absorber;
            Ec_abs = obj.params.chi_abs;
            Nc_abs = obj.params.Nc_abs;
            Ef(idx) = Ec_abs + obj.params.Vt * log(n(idx) ./ Nc_abs);
            
            % HTL区域
            idx = obj.params.idx_HTL;
            Ev_HTL = obj.params.chi_HTL - obj.params.Eg_HTL;
            Nv_HTL = obj.params.Nv_HTL;
            Ef(idx) = Ev_HTL - obj.params.Vt * log(p(idx) ./ Nv_HTL);
        end
        
        function setAppliedVoltage(obj, voltage)
            % 设置施加电压
            obj.params.setAppliedVoltage(voltage);
        end
        
        function converged = isConverged(obj, results, i, rel_tol, abs_tol)
            % 检查计算是否收敛
            
            % 对比当前时间步与前一步的结果
            phi_prev = results.phi(i-1,:);
            phi_curr = results.phi(i,:);
            
            n_prev = results.n(i-1,:);
            n_curr = results.n(i,:);
            
            p_prev = results.p(i-1,:);
            p_curr = results.p(i,:);
            
            % 计算相对变化
            phi_change = max(abs(phi_curr - phi_prev)) / (max(abs(phi_curr)) + abs_tol);
            n_change = max(abs(n_curr - n_prev)) / (max(abs(n_curr)) + abs_tol);
            p_change = max(abs(p_curr - p_prev)) / (max(abs(p_curr)) + abs_tol);
            
            % 当所有变量变化都小于容差时判定为收敛
            converged = (phi_change < rel_tol) && (n_change < rel_tol) && (p_change < rel_tol);
        end
    end
end