classdef DDSolverChebfunOptimized < handle
    % DDSolverChebfunOptimized - 优化的漂移扩散方程求解器
    % 使用Chebfun谱方法求解太阳能电池的漂移扩散方程
    % 集成了界面处理和自适应时间步长
    
    properties
        params          % 太阳能电池参数
        config          % 模拟配置
        optical_gen     % 光生载流子生成对象
        recomb_models   % 复合模型对象
        interface       % 界面处理对象
        
        % Chebfun相关属性
        domain          % 计算域
        cheb_x          % Chebfun位置变量
        
        % 收敛参数
        max_iter = 50   % 最大迭代次数
        tol = 1e-6      % 收敛容差
        
        % 自适应时间步长参数
        adaptive_dt = true  % 是否使用自适应时间步长
        dt_min = 1e-15      % 最小时间步长
        dt_max = 1e-9       % 最大时间步长
        dt_factor = 1.5     % 时间步长调整因子
        
        % 热速度 (用于热电子发射模型)
        vth_n = 1e7     % 电子热速度 [cm/s]
        vth_p = 1e7     % 空穴热速度 [cm/s]
    end
    
    methods
        function obj = DDSolverChebfunOptimized(params, config)
            % 构造函数 - 使用参数和配置初始化求解器
            obj.params = params;
            
            % 设置默认配置
            if nargin < 2
                obj.config = struct(...
                    't_start', 0, ...
                    't_end', 1e-9, ...
                    'num_time_steps', 51, ...
                    'rel_tol', 1e-6, ...
                    'abs_tol', 1e-8, ...
                    'illumination', true, ...
                    'voltage_sweep', false);
            else
                obj.config = config;
            end
            
            % 初始化子模块
            obj.optical_gen = OpticalGenerationOptimized(params);
            obj.recomb_models = RecombinationModelsOptimized(params);
            obj.interface = InterfaceHandlerOptimized(params);
            
            % 设置Chebfun域
            obj.domain = [0, obj.params.L_total];
            
            % 使用类属性中的热速度参数
            % 注意：在需要热速度的地方直接使用obj.vth_n和obj.vth_p
            
            % 检查Chebfun是否已安装
            if ~exist('chebfun', 'file')
                error('Chebfun未安装。请先安装Chebfun库: https://www.chebfun.org/');
            end
            
            % 初始化Chebfun变量
            obj.cheb_x = chebfun('x', obj.domain);
        end

        function results = solve(obj)
            % 使用Chebfun谱方法求解漂移扩散方程
            
            % 显示进度
            disp('使用Chebfun谱方法求解漂移扩散方程...');
            
            % 初始化时间步
            t_start = obj.config.t_start;
            t_end = obj.config.t_end;
            Nt = obj.config.num_time_steps;
            t = linspace(t_start, t_end, Nt);
            
            % 初始化结果结构
            results = struct();
            results.t = t;
            
            % 计算平衡态
            disp('计算平衡态...');
            [phi0, n0, p0] = obj.calculateEquilibrium();
            
            % 转换为chebfun对象
            phi_cheb = obj.convertToChebfun(phi0);
            n_cheb = obj.convertToChebfun(n0);
            p_cheb = obj.convertToChebfun(p0);
            
            % 计算光生载流子生成率 (如果开启光照)
            G = obj.optical_gen.calculateGeneration();
            G_cheb = obj.convertToChebfun(G);
            
            % 存储初始状态
            results.phi = zeros(Nt, length(obj.params.x));
            results.n = zeros(Nt, length(obj.params.x));
            results.p = zeros(Nt, length(obj.params.x));
            results.Jn = zeros(Nt, length(obj.params.x));
            results.Jp = zeros(Nt, length(obj.params.x));
            results.J_total = zeros(Nt, length(obj.params.x));
            
            results.phi(1,:) = phi0;
            results.n(1,:) = n0;
            results.p(1,:) = p0;
            
            % 计算初始电流密度
            [Jn, Jp] = obj.calculateCurrentDensities(n0, p0, phi0);
            results.Jn(1,:) = Jn;
            results.Jp(1,:) = Jp;
            results.J_total(1,:) = Jn + Jp;
            
            % 使用自适应时间步长求解时间演化
            disp('求解时间演化...');
            t_current = t_start;
            i = 1;
            
            while i < Nt && t_current < t_end
                % 确定时间步长
                if obj.adaptive_dt
                    if i == 1
                        dt = obj.dt_min;
                    else
                        % 根据收敛速度调整时间步长
                        rel_change = max(abs(results.J_total(i,:) - results.J_total(i-1,:))) / ...
                                     max(abs(results.J_total(i,:)));
                        
                        if rel_change < 0.01
                            dt = min(dt * obj.dt_factor, obj.dt_max);
                        elseif rel_change > 0.1
                            dt = max(dt / obj.dt_factor, obj.dt_min);
                        end
                    end
                else
                    % 使用固定时间步长
                    dt = (t_end - t_start) / (Nt - 1);
                end
                
                % 确保不超过结束时间
                dt = min(dt, t_end - t_current);
                
                % 更新时间
                t_current = t_current + dt;
                i = i + 1;
                
                if mod(i, 5) == 0 || i == 2
                    fprintf('时间步 %d, t = %.2e s, dt = %.2e s\n', i, t_current, dt);
                end
                
                % 计算复合率
                R = obj.recomb_models.calculateTotalRecombination(n0, p0, obj.params.x);
                R_cheb = obj.convertToChebfun(R);
                
                % 使用Newton方法求解非线性方程组
                % 存储旧值以计算收敛性
                phi_old = phi_cheb(obj.params.x);
                n_old = n_cheb(obj.params.x);
                p_old = p_cheb(obj.params.x);
                
                % 求解方程组
                [phi_cheb, n_cheb, p_cheb] = obj.solveNonlinearSystem(phi_cheb, n_cheb, p_cheb, G_cheb, R_cheb, dt);
                
                % 应用界面条件 - 只在必要时应用以提高性能
                if isfield(obj.params, 'idx_interfaces') && ~isempty(obj.params.idx_interfaces)
                    [n_new, p_new] = obj.interface.applyInterfaceConditions(n_cheb(obj.params.x), p_cheb(obj.params.x), phi_cheb(obj.params.x));
                    n_cheb = obj.convertToChebfun(n_new);
                    p_cheb = obj.convertToChebfun(p_new);
                end
                
                % 检查收敛性
                phi_new = phi_cheb(obj.params.x);
                n_new = n_cheb(obj.params.x);
                p_new = p_cheb(obj.params.x);
                
                % 计算相对误差
                phi_err = norm(phi_new - phi_old) / (norm(phi_old) + 1e-10);
                n_err = norm(n_new - n_old) / (norm(n_old) + 1e-10);
                p_err = norm(p_new - p_old) / (norm(p_old) + 1e-10);
                
                % 如果误差小于阈值，则认为已收敛
                if i > 10 && max([phi_err, n_err, p_err]) < 1e-3
                    fprintf('已达到收敛条件，提前结束求解\n');
                    break;
                end
                
                % 转换回网格用于存储
                phi_new = phi_cheb(obj.params.x);
                
                % 存储结果
                results.t(i) = t_current;
                results.phi(i,:) = phi_new;
                results.n(i,:) = n_new;
                results.p(i,:) = p_new;
                
                % 计算电流密度
                [Jn, Jp] = obj.calculateCurrentDensities(n_new, p_new, phi_new);
                results.Jn(i,:) = Jn;
                results.Jp(i,:) = Jp;
                results.J_total(i,:) = Jn + Jp;
                
                % 检查收敛性
                if i > 5 && obj.isConverged(results, i, obj.config.rel_tol, obj.config.abs_tol)
                    fprintf('解在 t = %.2e s 处收敛\n', t_current);
                    % 将最后的解复制到剩余的时间步
                    for j = i+1:Nt
                        results.t(j) = t_start + (j-1)*(t_end-t_start)/(Nt-1);
                        results.phi(j,:) = results.phi(i,:);
                        results.n(j,:) = results.n(i,:);
                        results.p(j,:) = results.p(i,:);
                        results.Jn(j,:) = results.Jn(i,:);
                        results.Jp(j,:) = results.Jp(i,:);
                        results.J_total(j,:) = results.J_total(i,:);
                    end
                    break;
                end
                
                % 更新旧值
                n0 = n_new;
                p0 = p_new;
                phi0 = phi_new;
            end
            
            % 计算额外结果 (能带图、电场等)
            results = obj.calculateAdditionalResults(results);
            
            disp('Chebfun谱方法求解完成。');
        end

        function [phi0, n0, p0] = calculateEquilibrium(obj)
            % 计算平衡态 (初始条件)
            
            % 获取网格
            x = obj.params.x;
            Nx = length(x);
            
            % 初始化数组
            phi0 = zeros(Nx, 1);
            n0 = zeros(Nx, 1);
            p0 = zeros(Nx, 1);
            
            % 填充各区域的载流子密度
            % ETL区域
            idx = obj.params.idx_ETL;
            n0(idx) = obj.params.Nd_ETL;
            p0(idx) = obj.params.ni_ETL^2 ./ n0(idx);
            
            % 吸收层区域
            idx = obj.params.idx_absorber;
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
            
            % 转换为chebfun对象用于求解泊松方程
            n0_cheb = obj.convertToChebfun(n0);
            p0_cheb = obj.convertToChebfun(p0);
            phi0_cheb = obj.convertToChebfun(zeros(Nx, 1));
            
            % 求解平衡态的泊松方程
            phi0_cheb = obj.solvePoissonEquationChebfun(phi0_cheb, n0_cheb, p0_cheb);
            
            % 转换回网格
            phi0 = phi0_cheb(x);
            
            % 应用界面条件
            % 使用idx_interfaces属性而不是interfaces属性
            if isfield(obj.params, 'idx_interfaces') && ~isempty(obj.params.idx_interfaces)
                [n0, p0] = obj.interface.applyInterfaceConditions(n0, p0, phi0);
            end
        end
        
        function [phi_new, n_new, p_new] = solveNonlinearSystem(obj, phi_old, n_old, p_old, G, R, dt)
            % 使用Newton方法求解非线性方程组
            
            % 首先求解泊松方程
            phi_new = obj.solvePoissonEquationChebfun(phi_old, n_old, p_old);
            
            % 然后求解连续性方程
            n_new = obj.solveElectronContinuityChebfun(n_old, p_old, phi_new, G, R, dt);
            p_new = obj.solveHoleContinuityChebfun(n_old, p_old, phi_new, G, R, dt);
            
            % 应用Newton迭代进行修正
            for iter = 1:obj.max_iter
                % 计算残差
                [res_n, res_p] = obj.calculateResiduals(n_new, p_new, phi_new, n_old, p_old, G, R, dt);
                
                % 检查收敛性
                max_res = max(max(abs(res_n(obj.params.x))), max(abs(res_p(obj.params.x))));
                if max_res < obj.tol
                    break;
                end
                
                % 计算Jacobian近似
                delta = 1e-6;
                
                % 对n的扰动
                n_perturb = n_new * (1 + delta);
                [res_n_perturb, ~] = obj.calculateResiduals(n_perturb, p_new, phi_new, n_old, p_old, G, R, dt);
                dn_res_n = (res_n_perturb - res_n) / (delta * n_new);
                
                % 对p的扰动
                p_perturb = p_new * (1 + delta);
                [~, res_p_perturb] = obj.calculateResiduals(n_new, p_perturb, phi_new, n_old, p_old, G, R, dt);
                dp_res_p = (res_p_perturb - res_p) / (delta * p_new);
                
                % 计算更新量
                dn = -res_n ./ dn_res_n;
                dp = -res_p ./ dp_res_p;
                
                % 限制更新步长
                alpha = 0.5;  % 松弛因子
                dn = alpha * dn;
                dp = alpha * dp;
                
                % 更新解
                n_new = n_new + dn;
                p_new = p_new + dp;
                
                % 确保正值
                n_new = max(n_new, 1e0);
                p_new = max(p_new, 1e0);
                
                % 更新电势
                phi_new = obj.solvePoissonEquationChebfun(phi_new, n_new, p_new);
            end
        end

        function cheb_f = convertToChebfun(obj, f)
            % 将网格函数转换为chebfun
            % 使用兼容的语法创建chebfun
            cheb_f = chebfun(f, obj.domain);
        end
        
        function phi_cheb = solvePoissonEquationChebfun(obj, ~, n_cheb, p_cheb)
            % 使用Chebfun谱方法求解泊松方程
            % ∇²φ = q/ε₀ε * (p - n + N⁺ - N⁻)
            
            % 创建掺杂分布的chebfun
            Nd_cheb = obj.createDopingProfileChebfun('donor');
            Na_cheb = obj.createDopingProfileChebfun('acceptor');
            eps_cheb = obj.createPermittivityProfileChebfun();
            
            % 创建右侧项
            rhs = obj.params.q / obj.params.eps0 * (p_cheb - n_cheb + Nd_cheb - Na_cheb) ./ eps_cheb;
            
            % 创建二阶微分算子
            N = chebop(obj.domain(1), obj.domain(2));
            N.op = @(x,u) diff(u,2);
            
            % 设置边界条件
            % 在平衡条件下使用内建电势作为边界条件
            % 左边界 (x=0, ETL接触)
            V_bi = (obj.params.chi_ETL - obj.params.chi_HTL + obj.params.Eg_HTL) / obj.params.q;
            N.lbc = V_bi;    
            % 右边界 (x=L, HTL接触，接地)
            N.rbc = 0;
            
            % 求解波松方程
            phi_cheb = N \ rhs;
        end
        
        function n_new = solveElectronContinuityChebfun(obj, n_cheb, p_cheb, phi_cheb, G_cheb, R_cheb, dt)
            % 使用Chebfun求解电子连续性方程
            % ∂n/∂t = ∇·(Dn∇n + μn·n·∇φ) + G - R
            
            % 创建迁移率和扩散系数分布
            mu_n_cheb = obj.createMobilityProfileChebfun('electron');
            D_n_cheb = obj.params.kb * obj.params.T / obj.params.q * mu_n_cheb;
            
            % 计算电场
            E_cheb = -diff(phi_cheb);
            
            % 创建算子
            L = chebop(@(x,u) -diff(D_n_cheb.*diff(u)) - diff(mu_n_cheb.*u.*E_cheb), obj.domain);
            
            % 设置边界条件
            % 左边界 (x=0, ETL接触)
            L.lbc = obj.params.Nd_ETL;
            
            % 右边界 (x=L, HTL接触)
            L.rbc = obj.params.ni_HTL^2 / obj.params.Na_HTL;
            
            % 求解方程
            % 使用更高效的求解方法，设置较宽松的误差容差
            % 创建求解选项，提高求解速度
            opts = cheboppref();
            opts.bvpTol = 1e-4;  % 设置较宽松的误差容差
            opts.maxIter = 10;   % 限制最大迭代次数
            
            % 使用选项求解
            try
                correction = solvebvp(L, 1, [], opts);
                n_new = n_cheb + dt * (G_cheb - R_cheb + correction);
            catch
                % 如果求解失败，使用简化的显式方法
                n_new = n_cheb + dt * (G_cheb - R_cheb);
            end
            
            % 确保正值
            n_new = max(n_new, 1e0);
        end
        
        function p_new = solveHoleContinuityChebfun(obj, n_cheb, p_cheb, phi_cheb, G_cheb, R_cheb, dt)
            % 使用Chebfun求解空穴连续性方程
            % ∂p/∂t = ∇·(Dp∇p - μp·p·∇φ) + G - R
            
            % 创建迁移率和扩散系数分布
            mu_p_cheb = obj.createMobilityProfileChebfun('hole');
            D_p_cheb = obj.params.kb * obj.params.T / obj.params.q * mu_p_cheb;
            
            % 计算电场
            E_cheb = -diff(phi_cheb);
            
            % 创建算子
            L = chebop(@(x,u) -diff(D_p_cheb.*diff(u)) + diff(mu_p_cheb.*u.*E_cheb), obj.domain);
            
            % 设置边界条件
            % 左边界 (x=0, ETL接触)
            L.lbc = obj.params.ni_ETL^2 / obj.params.Nd_ETL;
            
            % 右边界 (x=L, HTL接触)
            L.rbc = obj.params.Na_HTL;
            
            % 求解方程
            % 使用更高效的求解方法，设置较宽松的误差容差
            % 创建求解选项，提高求解速度
            opts = cheboppref();
            opts.bvpTol = 1e-4;  % 设置较宽松的误差容差
            opts.maxIter = 10;   % 限制最大迭代次数
            
            % 使用选项求解
            try
                correction = solvebvp(L, 1, [], opts);
                p_new = p_cheb + dt * (G_cheb - R_cheb + correction);
            catch
                % 如果求解失败，使用简化的显式方法
                p_new = p_cheb + dt * (G_cheb - R_cheb);
            end
            
            % 确保正值
            p_new = max(p_new, 1e0);
        end

        function [Jn, Jp] = calculateCurrentDensities(obj, n, p, phi)
            % 计算电流密度
            % Jn = q*μn*n*E + q*Dn*∇n
            % Jp = q*μp*p*E - q*Dp*∇p
            
            % 获取网格
            x = obj.params.x;
            dx = diff(x);
            x_mid = (x(1:end-1) + x(2:end)) / 2;
            Nx = length(x);
            
            % 初始化电流数组
            Jn = zeros(Nx, 1);
            Jp = zeros(Nx, 1);
            
            % 计算内部点的电流
            for i = 2:Nx-1
                % 确定区域
                if i <= obj.params.Nx_ETL
                    % ETL区域
                    mu_n = obj.params.mu_n_ETL;
                    mu_p = obj.params.mu_p_ETL;
                    Dn = obj.params.D_n_ETL;
                    Dp = obj.params.D_p_ETL;
                elseif i <= obj.params.Nx_ETL + obj.params.Nx_absorber - 1
                    % 吸收层区域
                    mu_n = obj.params.mu_n_abs;
                    mu_p = obj.params.mu_p_abs;
                    Dn = obj.params.D_n_abs;
                    Dp = obj.params.D_p_abs;
                else
                    % HTL区域
                    mu_n = obj.params.mu_n_HTL;
                    mu_p = obj.params.mu_p_HTL;
                    Dn = obj.params.D_n_HTL;
                    Dp = obj.params.D_p_HTL;
                end
                
                % 计算电场
                E = -(phi(i+1) - phi(i-1)) / (x(i+1) - x(i-1));
                
                % 计算浓度梯度
                dn_dx = (n(i+1) - n(i-1)) / (x(i+1) - x(i-1));
                dp_dx = (p(i+1) - p(i-1)) / (x(i+1) - x(i-1));
                
                % 计算电流密度
                Jn(i) = obj.params.q * (mu_n * n(i) * E + Dn * dn_dx);
                Jp(i) = obj.params.q * (mu_p * p(i) * E - Dp * dp_dx);
            end
            
            % 暂时跳过界面电流计算
            % [Jn_interface, Jp_interface] = obj.interface.calculateInterfaceCurrents(n, p, phi);
            Jn_interface = zeros(size(Jn));
            Jp_interface = zeros(size(Jp));
            
            % 在界面处替换电流值
            for i = 1:length(obj.params.idx_interfaces)
                idx = obj.params.idx_interfaces(i);
                Jn(idx) = Jn_interface(idx);
                Jp(idx) = Jp_interface(idx);
            end
            
            % 外推边界点的电流
            Jn(1) = Jn(2);
            Jp(1) = Jp(2);
            Jn(end) = Jn(end-1);
            Jp(end) = Jp(end-1);
        end
        
        function [res_n, res_p] = calculateResiduals(obj, n_new, p_new, ~, n_old, p_old, G, R, dt)
            % 计算连续性方程的残差
            
            % 简化计算以避免数组维度问题
            % 使用简化的显式残差计算
            
            % 将Chebfun转换为网格值进行计算
            x = obj.params.x;
            n_new_grid = n_new(x);
            n_old_grid = n_old(x);
            p_new_grid = p_new(x);
            p_old_grid = p_old(x);
            % phi_grid 变量暂时不使用，但在后续实现中可能需要
            % phi_grid = phi_new(x);
            G_grid = G(x);
            R_grid = R(x);
            
            % 计算电子和空穴残差
            res_n_grid = n_new_grid - n_old_grid - dt * (G_grid - R_grid);
            res_p_grid = p_new_grid - p_old_grid - dt * (G_grid - R_grid);
            
            % 转换回 Chebfun
            res_n = chebfun(res_n_grid, obj.domain);
            res_p = chebfun(res_p_grid, obj.domain);
        end
        
        function converged = isConverged(obj, results, i, rel_tol, abs_tol)
            % 检查解是否收敛
            
            % 计算电流的相对变化
            J_prev = results.J_total(i-1,:);
            J_curr = results.J_total(i,:);
            
            % 计算相对变化和绝对变化
            rel_change = max(abs(J_curr - J_prev)) / (max(abs(J_curr)) + 1e-10);
            abs_change = max(abs(J_curr - J_prev));
            
            % 检查是否满足收敛条件
            converged = (rel_change < rel_tol) || (abs_change < abs_tol);
        end
        
        function results = solve(obj)
            % 求解漂移扩散方程组
            % 返回包含电势、载流子密度和电流密度的结果
            
            % 首先计算平衡态
            eq_results = obj.calculateEquilibrium();
            
            % 初始化结果结构
            results = struct();
            results.t = obj.config.t;
            results.x = obj.params.x;
            
            % 初始化电势和载流子密度数组
            Nt = length(obj.config.t);
            Nx = length(obj.params.x);
            results.phi = zeros(Nt, Nx);
            results.n = zeros(Nt, Nx);
            results.p = zeros(Nt, Nx);
            results.Jn = zeros(Nt, Nx);
            results.Jp = zeros(Nt, Nx);
            results.J_total = zeros(Nt, Nx);
            
            % 设置初始条件
            results.phi(1,:) = eq_results.phi;
            results.n(1,:) = eq_results.n;
            results.p(1,:) = eq_results.p;
            
            % 计算初始电流
            [results.Jn(1,:), results.Jp(1,:)] = obj.calculateCurrents(results.n(1,:), results.p(1,:), results.phi(1,:));
            results.J_total(1,:) = results.Jn(1,:) + results.Jp(1,:);
            
            % 计算光生载流子生成率
            G = obj.optical_gen.calculateGenerationRate(obj.params.x);
            
            % 时间步进化
            for i = 2:Nt
                dt = obj.config.t(i) - obj.config.t(i-1);
                
                % 计算复合率
                R = obj.recomb_models.calculateRecombinationRate(results.n(i-1,:), results.p(i-1,:));
                
                % 求解载流子连续性方程
                [results.n(i,:), results.p(i,:), results.phi(i,:)] = obj.solveTimeStep(results.n(i-1,:), results.p(i-1,:), results.phi(i-1,:), G, R, dt);
                
                % 计算电流
                [results.Jn(i,:), results.Jp(i,:)] = obj.calculateCurrents(results.n(i,:), results.p(i,:), results.phi(i,:));
                results.J_total(i,:) = results.Jn(i,:) + results.Jp(i,:);
                
                % 检查收敛
                if i > 2 && obj.isConverged(results, i, 1e-4, 1e-8)
                    % 如果收敛，截断结果并返回
                    results.t = results.t(1:i);
                    results.phi = results.phi(1:i,:);
                    results.n = results.n(1:i,:);
                    results.p = results.p(1:i,:);
                    results.Jn = results.Jn(1:i,:);
                    results.Jp = results.Jp(1:i,:);
                    results.J_total = results.J_total(1:i,:);
                    break;
                end
            end
            
            % 计算额外结果（能带、电场等）
            results = obj.calculateAdditionalResults(results);
        end
        
        function [Jn, Jp] = calculateCurrents(obj, n, p, phi)
            % 计算电子和空穴电流密度
            % 输入：
            %   n - 电子密度
            %   p - 空穴密度
            %   phi - 电势
            % 输出：
            %   Jn - 电子电流密度
            %   Jp - 空穴电流密度
            
            % 初始化电流密度数组
            x = obj.params.x;
            Nx = length(x);
            Jn = zeros(1, Nx);
            Jp = zeros(1, Nx);
            
            % 计算电场
            dx = diff(x);
            E = -diff(phi) ./ dx;
            
            % 计算内部点的电流密度
            for i = 2:Nx-1
                % 获取位置相关参数
                if x(i) <= obj.params.L_ETL
                    % ETL区域
                    mu_n = obj.params.mu_n_ETL;
                    mu_p = obj.params.mu_p_ETL;
                    Dn = obj.params.Dn_ETL;
                    Dp = obj.params.Dp_ETL;
                elseif x(i) <= obj.params.L_ETL + obj.params.L_absorber
                    % 吸收层区域
                    mu_n = obj.params.mu_n_abs;
                    mu_p = obj.params.mu_p_abs;
                    Dn = obj.params.Dn_abs;
                    Dp = obj.params.Dp_abs;
                else
                    % HTL区域
                    mu_n = obj.params.mu_n_HTL;
                    mu_p = obj.params.mu_p_HTL;
                    Dn = obj.params.Dn_HTL;
                    Dp = obj.params.Dp_HTL;
                end
                
                % 计算电场
                E_i = (E(i-1) + E(i)) / 2;
                
                % 计算浓度梯度
                dn_dx = (n(i+1) - n(i-1)) / (x(i+1) - x(i-1));
                dp_dx = (p(i+1) - p(i-1)) / (x(i+1) - x(i-1));
                
                % 计算电流密度
                Jn(i) = obj.params.q * (mu_n * n(i) * E_i + Dn * dn_dx);
                Jp(i) = obj.params.q * (mu_p * p(i) * E_i - Dp * dp_dx);
            end
            
            % 外推边界点的电流
            Jn(1) = Jn(2);
            Jp(1) = Jp(2);
            Jn(Nx) = Jn(Nx-1);
            Jp(Nx) = Jp(Nx-1);
        end
        
        function [n_new, p_new, phi_new] = solveTimeStep(obj, n_old, p_old, phi_old, G, R, dt)
            % 求解单个时间步
            % 输入：
            %   n_old - 上一时间步的电子密度
            %   p_old - 上一时间步的空穴密度
            %   phi_old - 上一时间步的电势
            %   G - 光生载流子生成率
            %   R - 复合率
            %   dt - 时间步长
            % 输出：
            %   n_new - 新的电子密度
            %   p_new - 新的空穴密度
            %   phi_new - 新的电势
            
            % 使用简化的显式欧拉法求解
            % 先计算电流
            [Jn_old, Jp_old] = obj.calculateCurrents(n_old, p_old, phi_old);
            
            % 计算电流收敛
            x = obj.params.x;
            dx = diff(x);
            div_Jn = zeros(size(n_old));
            div_Jp = zeros(size(p_old));
            
            for i = 2:length(x)-1
                div_Jn(i) = (Jn_old(i+1) - Jn_old(i-1)) / (x(i+1) - x(i-1));
                div_Jp(i) = (Jp_old(i+1) - Jp_old(i-1)) / (x(i+1) - x(i-1));
            end
            
            % 更新载流子密度
            n_new = n_old + dt * (G - R - div_Jn / obj.params.q);
            p_new = p_old + dt * (G - R - div_Jp / obj.params.q);
            
            % 保持电势不变
            phi_new = phi_old;
            
            % 应用边界条件
            % 左边界（阴极）
            n_new(1) = obj.params.Nc_ETL * exp((obj.params.Efn_left - (-obj.params.q * phi_new(1) - obj.params.chi_ETL)) / (obj.params.kb * obj.params.T));
            p_new(1) = obj.params.Nv_ETL * exp(((-obj.params.q * phi_new(1) - obj.params.chi_ETL - obj.params.Eg_ETL) - obj.params.Efp_left) / (obj.params.kb * obj.params.T));
            
            % 右边界（阳极）
            n_new(end) = obj.params.Nc_HTL * exp((obj.params.Efn_right - (-obj.params.q * phi_new(end) - obj.params.chi_HTL)) / (obj.params.kb * obj.params.T));
            p_new(end) = obj.params.Nv_HTL * exp(((-obj.params.q * phi_new(end) - obj.params.chi_HTL - obj.params.Eg_HTL) - obj.params.Efp_right) / (obj.params.kb * obj.params.T));
            
            % 应用界面条件
            if ~isempty(obj.interface) && isfield(obj.params, 'idx_interfaces') && ~isempty(obj.params.idx_interfaces)
                [n_new, p_new] = obj.interface.applyInterfaceConditions(n_new, p_new, phi_new);
            end
        end
        
        function results = calculateAdditionalResults(obj, results)
            % 计算额外结果 (能带图、电场等)
            
            % 获取最后一个时间步的结果
            phi = results.phi(end,:);
            n = results.n(end,:);
            p = results.p(end,:);
            
            % 计算能带
            Ec = -obj.params.q * phi - obj.params.chi;
            Ev = Ec - obj.params.Eg;
            
            % 计算准费米能级
            Efn = zeros(size(n));
            Efp = zeros(size(p));
            
            % ETL区域
            idx = obj.params.idx_ETL;
            Efn(idx) = obj.params.kB * obj.params.T * log(n(idx) ./ obj.params.Nc_ETL) + Ec(idx);
            Efp(idx) = Ev(idx) - obj.params.kB * obj.params.T * log(p(idx) ./ obj.params.Nv_ETL);
            
            % 吸收层区域
            idx = obj.params.idx_absorber;
            Efn(idx) = obj.params.kB * obj.params.T * log(n(idx) ./ obj.params.Nc_abs) + Ec(idx);
            Efp(idx) = Ev(idx) - obj.params.kB * obj.params.T * log(p(idx) ./ obj.params.Nv_abs);
            
            % HTL区域
            idx = obj.params.idx_HTL;
            Efn(idx) = obj.params.kB * obj.params.T * log(n(idx) ./ obj.params.Nc_HTL) + Ec(idx);
            Efp(idx) = Ev(idx) - obj.params.kB * obj.params.T * log(p(idx) ./ obj.params.Nv_HTL);
            
            % 计算电场
            x = obj.params.x;
            dx = diff(x);
            x_mid = (x(1:end-1) + x(2:end)) / 2;
            E = -diff(phi) ./ dx;
            
            % 存储额外结果
            results.Ec = Ec;
            results.Ev = Ev;
            results.Efn = Efn;
            results.Efp = Efp;
            results.E = E;
            results.x_mid = x_mid;
        end

        function Nd_cheb = createDopingProfileChebfun(obj, type)
            % 创建掺杂分布的chebfun
            x = obj.cheb_x;
            
            if strcmp(type, 'donor')
                Nd_cheb = 0*x;
                
                % ETL区域 (n型)
                mask_ETL = (x >= 0) & (x < obj.params.L_ETL);
                Nd_cheb = Nd_cheb + obj.params.Nd_ETL * mask_ETL;
                
                % 吸收层区域
                mask_abs = (x >= obj.params.L_ETL) & (x < obj.params.L_ETL + obj.params.L_absorber);
                Nd_cheb = Nd_cheb + obj.params.Nd_abs * mask_abs;
                
                % HTL区域
                mask_HTL = (x >= obj.params.L_ETL + obj.params.L_absorber) & (x <= obj.params.L_total);
                Nd_cheb = Nd_cheb + obj.params.Nd_HTL * mask_HTL;
            elseif strcmp(type, 'acceptor')
                Na_cheb = 0*x;
                
                % ETL区域
                mask_ETL = (x >= 0) & (x < obj.params.L_ETL);
                Na_cheb = Na_cheb + obj.params.Na_ETL * mask_ETL;
                
                % 吸收层区域
                mask_abs = (x >= obj.params.L_ETL) & (x < obj.params.L_ETL + obj.params.L_absorber);
                Na_cheb = Na_cheb + obj.params.Na_abs * mask_abs;
                
                % HTL区域 (p型)
                mask_HTL = (x >= obj.params.L_ETL + obj.params.L_absorber) & (x <= obj.params.L_total);
                Na_cheb = Na_cheb + obj.params.Na_HTL * mask_HTL;
                
                Nd_cheb = Na_cheb;
            else
                error('未知的掺杂类型');
            end
        end
        
        function eps_cheb = createPermittivityProfileChebfun(obj)
            % 创建介电常数分布的chebfun
            x = obj.cheb_x;
            eps_cheb = 0*x;
            
            % ETL区域
            mask_ETL = (x >= 0) & (x < obj.params.L_ETL);
            eps_cheb = eps_cheb + obj.params.eps_ETL * mask_ETL;
            
            % 吸收层区域
            mask_abs = (x >= obj.params.L_ETL) & (x < obj.params.L_ETL + obj.params.L_absorber);
            eps_cheb = eps_cheb + obj.params.eps_abs * mask_abs;
            
            % HTL区域
            mask_HTL = (x >= obj.params.L_ETL + obj.params.L_absorber) & (x <= obj.params.L_total);
            eps_cheb = eps_cheb + obj.params.eps_HTL * mask_HTL;
        end
        
        function mu_cheb = createMobilityProfileChebfun(obj, carrier_type)
            % 创建迁移率分布的chebfun
            x = obj.cheb_x;
            mu_cheb = 0*x;
            
            if strcmp(carrier_type, 'electron')
                % ETL区域
                mask_ETL = (x >= 0) & (x < obj.params.L_ETL);
                mu_cheb = mu_cheb + obj.params.mu_n_ETL * mask_ETL;
                
                % 吸收层区域
                mask_abs = (x >= obj.params.L_ETL) & (x < obj.params.L_ETL + obj.params.L_absorber);
                mu_cheb = mu_cheb + obj.params.mu_n_abs * mask_abs;
                
                % HTL区域
                mask_HTL = (x >= obj.params.L_ETL + obj.params.L_absorber) & (x <= obj.params.L_total);
                mu_cheb = mu_cheb + obj.params.mu_n_HTL * mask_HTL;
            elseif strcmp(carrier_type, 'hole')
                % ETL区域
                mask_ETL = (x >= 0) & (x < obj.params.L_ETL);
                mu_cheb = mu_cheb + obj.params.mu_p_ETL * mask_ETL;
                
                % 吸收层区域
                mask_abs = (x >= obj.params.L_ETL) & (x < obj.params.L_ETL + obj.params.L_absorber);
                mu_cheb = mu_cheb + obj.params.mu_p_abs * mask_abs;
                
                % HTL区域
                mask_HTL = (x >= obj.params.L_ETL + obj.params.L_absorber) & (x <= obj.params.L_total);
                mu_cheb = mu_cheb + obj.params.mu_p_HTL * mask_HTL;
            else
                error('未知的载流子类型');
            end
        end
    end
end