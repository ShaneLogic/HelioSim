classdef VisualizerOptimized < handle
    % VisualizerOptimized - 优化的结果可视化工具
    % 提供能带图、载流子密度、电场分布和J-V曲线的可视化
    
    properties
        params      % 太阳能电池参数
        results     % 模拟结果
        jv_results  % J-V分析结果
        fig_handles % 图形句柄
    end
    
    methods
        function obj = VisualizerOptimized(params)
            % 构造函数 - 初始化可视化工具
            obj.params = params;
            obj.fig_handles = struct();
        end
        
        function setResults(obj, results)
            % 设置模拟结果
            obj.results = results;
        end
        
        function setJVResults(obj, jv_results)
            % 设置J-V分析结果
            obj.jv_results = jv_results;
        end
        
        function plotBandDiagram(obj)
            % 绘制能带图
            
            % 检查结果是否存在
            if isempty(obj.results)
                error('没有可用的模拟结果。请先运行模拟。');
            end
            
            % 获取数据
            x = obj.results.x;  % 获取x坐标数据
            
            % 检查数据结构并根据需要调整维度
            if isfield(obj.results, 't') && length(obj.results.t) > 1
                % 多时间步结果，获取最后一个时间步
                phi = obj.results.phi(end, :);
                n = obj.results.n(end, :);
                p = obj.results.p(end, :);
            else
                % 单个时间步结果
                phi = obj.results.phi;
                n = obj.results.n;
                p = obj.results.p;
            end
            
            % 计算能带
            % 创建电子亲和能分布
            chi = zeros(size(x));
            Eg = zeros(size(x));
            Nc = zeros(size(x));
            Nv = zeros(size(x));
            
            % 根据位置设置不同区域的参数
            for i = 1:length(x)
                if x(i) <= obj.params.L_ETL
                    % ETL区域
                    chi(i) = obj.params.chi_ETL;
                    Eg(i) = obj.params.Eg_ETL;
                    Nc(i) = obj.params.Nc_ETL;
                    Nv(i) = obj.params.Nv_ETL;
                elseif x(i) <= obj.params.L_ETL + obj.params.L_absorber
                    % 吸收层区域
                    chi(i) = obj.params.chi_abs;
                    Eg(i) = obj.params.Eg_abs;
                    Nc(i) = obj.params.Nc_abs;
                    Nv(i) = obj.params.Nv_abs;
                else
                    % HTL区域
                    chi(i) = obj.params.chi_HTL;
                    Eg(i) = obj.params.Eg_HTL;
                    Nc(i) = obj.params.Nc_HTL;
                    Nv(i) = obj.params.Nv_HTL;
                end
            end
            
            % 计算能带
            Ec = -obj.params.q * phi - chi;
            Ev = Ec - Eg;
            Efn = Ec + obj.params.kb * obj.params.T * log(n ./ Nc);
            Efp = Ev - obj.params.kb * obj.params.T * log(p ./ Nv);
            
            % 创建图形
            obj.fig_handles.band = figure('Name', '能带图', 'Position', [100, 100, 800, 600]);
            
            % 绘制能带
            plot(obj.params.x*1e7, Ec, 'b-', 'LineWidth', 2);
            hold on;
            plot(obj.params.x*1e7, Ev, 'r-', 'LineWidth', 2);
            plot(obj.params.x*1e7, Efn, 'b--', 'LineWidth', 1.5);
            plot(obj.params.x*1e7, Efp, 'r--', 'LineWidth', 1.5);
            
            % 绘制界面
            % 使用idx_interfaces而不interfaces
            if isfield(obj.params, 'idx_interfaces') && ~isempty(obj.params.idx_interfaces)
                interface_positions = [obj.params.L_ETL, obj.params.L_ETL + obj.params.L_absorber];
                for i = 1:length(interface_positions)
                    x_int = interface_positions(i);
                    plot([x_int, x_int]*1e7, [min(Ev)-0.5, max(Ec)+0.5], 'k--');
                end
            end
            
            % 添加层标签
            layers = {'ETL', '吸收层', 'HTL'};
            x_centers = zeros(1, 3);
            x_centers(1) = mean(obj.params.x(obj.params.idx_ETL))*1e7;
            x_centers(2) = mean(obj.params.x(obj.params.idx_absorber))*1e7;
            x_centers(3) = mean(obj.params.x(obj.params.idx_HTL))*1e7;
            
            for i = 1:length(layers)
                text(x_centers(i), max(Ec)+0.3, layers{i}, ...
                     'HorizontalAlignment', 'center', 'FontWeight', 'bold');
            end
            
            % 设置图形属性
            xlabel('位置 (nm)');
            ylabel('能量 (eV)');
            title('能带图');
            legend('导带', '价带', '电子准费米能级', '空穴准费米能级', 'Location', 'best');
            grid on;
            
            % 设置坐标轴范围
            xlim([min(obj.params.x), max(obj.params.x)]*1e7);
            ylim([min(Ev)-0.5, max(Ec)+0.5]);
            
            % 反转Y轴
            set(gca, 'YDir', 'reverse');
        end
        
        function plotCarrierDensities(obj)
            % 绘制载流子密度分布
            
            % 检查结果是否存在
            if isempty(obj.results)
                error('没有可用的模拟结果。请先运行模拟。');
            end
            
            % 获取数据
            % 检查数据结构并根据需要调整维度
            if isfield(obj.results, 't') && length(obj.results.t) > 1
                % 多时间步结果，获取最后一个时间步
                n = obj.results.n(end, :);
                p = obj.results.p(end, :);
            else
                % 单个时间步结果
                n = obj.results.n;
                p = obj.results.p;
            end
            
            % 创建图形
            obj.fig_handles.carriers = figure('Name', '载流子密度', 'Position', [100, 100, 800, 600]);
            
            % 绘制载流子密度
            semilogy(obj.params.x*1e7, n, 'b-', 'LineWidth', 2);
            hold on;
            semilogy(obj.params.x*1e7, p, 'r-', 'LineWidth', 2);
            
            % 添加界面线
            % 使用固定的界面位置
            interface_positions = [obj.params.L_ETL, obj.params.L_ETL + obj.params.L_absorber];
            for i = 1:length(interface_positions)
                x_int = interface_positions(i);
                plot([x_int, x_int]*1e7, [1e0, 1e20], 'k--');
            end
            
            % 添加层标签
            layers = {'ETL', '吸收层', 'HTL'};
            x_centers = zeros(1, 3);
            x_centers(1) = mean(obj.params.x(obj.params.idx_ETL))*1e7;
            x_centers(2) = mean(obj.params.x(obj.params.idx_absorber))*1e7;
            x_centers(3) = mean(obj.params.x(obj.params.idx_HTL))*1e7;
            
            for i = 1:length(layers)
                text(x_centers(i), 1e19, layers{i}, ...
                     'HorizontalAlignment', 'center', 'FontWeight', 'bold');
            end
            
            % 设置图形属性
            xlabel('位置 (nm)');
            ylabel('载流子密度 (cm^{-3})');
            title('载流子密度分布');
            legend('电子密度', '空穴密度', 'Location', 'best');
            grid on;
            
            % 设置坐标轴范围
            xlim([min(obj.params.x), max(obj.params.x)]*1e7);
            ylim([1e0, 1e20]);
        end
        
        function plotElectricField(obj)
            % 绘制电场分布
            
            % 检查结果是否存在
            if isempty(obj.results)
                error('没有可用的模拟结果。请先运行模拟。');
            end
            
            % 获取数据
            % 检查数据结构并根据需要调整维度
            if isfield(obj.results, 't') && length(obj.results.t) > 1
                % 多时间步结果，获取最后一个时间步
                phi = obj.results.phi(end, :);
            else
                % 单个时间步结果
                phi = obj.results.phi;
            end
            
            % 计算电场
            x = obj.params.x;
            dx = diff(x);
            x_mid = (x(1:end-1) + x(2:end)) / 2;
            E = -diff(phi) ./ dx;  % V/cm
            
            % 创建图形
            obj.fig_handles.field = figure('Name', '电场分布', 'Position', [100, 100, 800, 600]);
            
            % 绘制电场
            plot(x_mid*1e7, E, 'k-', 'LineWidth', 2);
            
            % 添加界面线
            % 使用固定的界面位置
            interface_positions = [obj.params.L_ETL, obj.params.L_ETL + obj.params.L_absorber];
            for i = 1:length(interface_positions)
                x_int = interface_positions(i);
                y_lim = get(gca, 'YLim');
                plot([x_int, x_int]*1e7, y_lim, 'k--');
            end
            
            % 添加层标签
            layers = {'ETL', '吸收层', 'HTL'};
            x_centers = zeros(1, 3);
            x_centers(1) = mean(obj.params.x(obj.params.idx_ETL))*1e7;
            x_centers(2) = mean(obj.params.x(obj.params.idx_absorber))*1e7;
            x_centers(3) = mean(obj.params.x(obj.params.idx_HTL))*1e7;
            
            for i = 1:length(layers)
                text(x_centers(i), max(E)*0.9, layers{i}, ...
                     'HorizontalAlignment', 'center', 'FontWeight', 'bold');
            end
            
            % 设置图形属性
            xlabel('位置 (nm)');
            ylabel('电场 (V/cm)');
            title('电场分布');
            grid on;
            
            % 设置坐标轴范围
            xlim([min(obj.params.x), max(obj.params.x)]*1e7);
        end
        
        function plotCurrentDensities(obj)
            % 绘制电流密度分布
            
            % 检查结果是否存在
            if isempty(obj.results)
                error('没有可用的模拟结果。请先运行模拟。');
            end
            
            % 获取最后一个时间步的结果
            if isfield(obj.results, 't') && length(obj.results.t) > 1
                % 多时间步结果，获取最后一个时间步
                Jn = obj.results.Jn(end, :);
                Jp = obj.results.Jp(end, :);
            else
                % 单个时间步结果
                Jn = obj.results.Jn;
                Jp = obj.results.Jp;
            end
            Jtot = Jn + Jp;
            
            % 创建图形
            obj.fig_handles.current = figure('Name', '电流密度', 'Position', [100, 100, 800, 600]);
            
            % 绘制电流密度
            plot(obj.params.x*1e7, Jn, 'b-', 'LineWidth', 2);
            hold on;
            plot(obj.params.x*1e7, Jp, 'r-', 'LineWidth', 2);
            plot(obj.params.x*1e7, Jtot, 'k-', 'LineWidth', 2);
            
            % 添加界面线
            % 使用固定的界面位置
            interface_positions = [obj.params.L_ETL, obj.params.L_ETL + obj.params.L_absorber];
            for i = 1:length(interface_positions)
                x_int = interface_positions(i);
                y_lim = get(gca, 'YLim');
                plot([x_int, x_int]*1e7, y_lim, 'k--');
            end
            
            % 添加层标签
            layers = {'ETL', '吸收层', 'HTL'};
            x_centers = zeros(1, 3);
            x_centers(1) = mean(obj.params.x(obj.params.idx_ETL))*1e7;
            x_centers(2) = mean(obj.params.x(obj.params.idx_absorber))*1e7;
            x_centers(3) = mean(obj.params.x(obj.params.idx_HTL))*1e7;
            
            for i = 1:length(layers)
                text(x_centers(i), max(Jtot)*0.9, layers{i}, ...
                     'HorizontalAlignment', 'center', 'FontWeight', 'bold');
            end
            
            % 设置图形属性
            xlabel('位置 (nm)');
            ylabel('电流密度 (mA/cm^2)');
            title('电流密度分布');
            legend('电子电流', '空穴电流', '总电流', 'Location', 'best');
            grid on;
            
            % 设置坐标轴范围
            xlim([min(obj.params.x), max(obj.params.x)]*1e7);
        end
        
        function plotJVCurve(obj)
            % 绘制J-V特性曲线
            
            % 检查J-V结果是否存在
            if isempty(obj.jv_results)
                error('没有可用的J-V分析结果。请先运行J-V扫描。');
            end
            
            % 创建图形
            obj.fig_handles.jv = figure('Name', 'J-V特性曲线', 'Position', [100, 100, 800, 600]);
            
            % 绘制J-V曲线
            plot(obj.jv_results.V, obj.jv_results.J, 'k-', 'LineWidth', 2);
            
            % 标记关键点
            hold on;
            plot(obj.jv_results.Voc, 0, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
            plot(0, obj.jv_results.Jsc, 'bo', 'MarkerSize', 10, 'LineWidth', 2);
            plot(obj.jv_results.V_mpp, obj.jv_results.J_mpp, 'go', 'MarkerSize', 10, 'LineWidth', 2);
            
            % 设置图形属性
            xlabel('电压 (V)');
            ylabel('电流密度 (mA/cm^2)');
            title('J-V特性曲线');
            grid on;
            
            % 添加关键参数文本
            text_str = sprintf('Voc = %.3f V\nJsc = %.2f mA/cm^2\nFF = %.2f %%\nPCE = %.2f %%', ...
                          obj.jv_results.Voc, obj.jv_results.Jsc, ...
                          obj.jv_results.FF*100, obj.jv_results.PCE);
            text(0.05, 0.95, text_str, 'Units', 'normalized', ...
                 'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');
        end
        
        function plotHysteresisJVCurve(obj)
            % 绘制带迟滞的J-V特性曲线
            
            % 检查J-V结果是否存在
            if ~isfield(obj.jv_results, 'forward') || ~isfield(obj.jv_results, 'reverse')
                error('没有可用的迟滞J-V分析结果。请先运行正向和反向J-V扫描。');
            end
            
            % 创建图形
            obj.fig_handles.jv_hysteresis = figure('Name', '迟滞J-V特性曲线', 'Position', [100, 100, 800, 600]);
            
            % 绘制正向扫描J-V曲线
            plot(obj.jv_results.forward.V, obj.jv_results.forward.J, 'b-', 'LineWidth', 2);
            hold on;
            
            % 绘制反向扫描J-V曲线
            plot(obj.jv_results.reverse.V, obj.jv_results.reverse.J, 'r-', 'LineWidth', 2);
            
            % 设置图形属性
            xlabel('电压 (V)');
            ylabel('电流密度 (mA/cm^2)');
            title('迟滞J-V特性曲线');
            legend('正向扫描', '反向扫描', 'Location', 'best');
            grid on;
            
            % 添加关键参数文本
            text_str = sprintf('正向扫描:\nVoc = %.3f V\nJsc = %.2f mA/cm^2\nFF = %.2f %%\nPCE = %.2f %%\n\n反向扫描:\nVoc = %.3f V\nJsc = %.2f mA/cm^2\nFF = %.2f %%\nPCE = %.2f %%', ...
                          obj.jv_results.forward.Voc, obj.jv_results.forward.Jsc, ...
                          obj.jv_results.forward.FF*100, obj.jv_results.forward.PCE, ...
                          obj.jv_results.reverse.Voc, obj.jv_results.reverse.Jsc, ...
                          obj.jv_results.reverse.FF*100, obj.jv_results.reverse.PCE);
            text(0.05, 0.95, text_str, 'Units', 'normalized', ...
                 'VerticalAlignment', 'top', 'FontSize', 10, 'FontWeight', 'bold');
        end
    end
end
