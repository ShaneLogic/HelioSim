classdef JVAnalyzerOptimized < handle
    % JVAnalyzerOptimized - 优化的J-V特性曲线分析器
    % 计算太阳能电池的性能参数，包括Voc、Jsc、FF和效率
    
    properties
        params          % 太阳能电池参数
        solver          % 漂移扩散求解器
        AM15_power = 100; % AM1.5光谱功率密度 [mW/cm^2]
    end
    
    methods
        function obj = JVAnalyzerOptimized(params, solver)
            % 构造函数 - 初始化分析器
            obj.params = params;
            obj.solver = solver;
        end
        
        function jv_results = generateJVCurve(obj, V_start, V_end, num_points)
            % 生成J-V特性曲线
            % V_start: 起始电压 [V]
            % V_end: 终止电压 [V]
            % num_points: 电压点数
            
            if nargin < 4
                num_points = 30; % 默认点数
            end
            
            if nargin < 3
                V_start = -0.2; % 默认起始电压
                V_end = 1.2;    % 默认终止电压
            end
            
            % 确保光照开启
            illumination_state = obj.params.illumination;
            obj.params.setIllumination(true);
            
            % 创建电压数组
            voltages = linspace(V_start, V_end, num_points);
            
            % 初始化电流密度数组
            currents = zeros(size(voltages));
            
            % 计算每个电压点的电流密度
            disp('生成J-V曲线...');
            for i = 1:length(voltages)
                fprintf('计算电压点 %d/%d: %.3f V\n', i, length(voltages), voltages(i));
                
                % 设置施加电压
                obj.params.setAppliedVoltage(voltages(i));
                
                % 求解漂移扩散方程
                results = obj.solver.solve();
                
                % 提取电流密度 (mA/cm^2)
                % 使用最后一个时间步的总电流
                J_total = results.J_total(end, :);
                
                % 使用接触处的电流密度
                currents(i) = J_total(1) * 1e3; % 转换为mA/cm^2
            end
            
            % 恢复原始光照状态
            obj.params.setIllumination(illumination_state);
            
            % 计算性能参数
            jv_results = obj.calculatePerformanceParameters(voltages, currents);
            
            % 存储J-V数据
            jv_results.voltages = voltages;
            jv_results.currents = currents;
        end
        
        function params = calculatePerformanceParameters(obj, voltages, currents)
            % 计算太阳能电池性能参数
            % 包括Voc、Jsc、FF和效率
            
            % 初始化结果结构
            params = struct();
            
            % 计算短路电流 (mA/cm^2)
            % 找到最接近V=0的点
            [~, idx_sc] = min(abs(voltages));
            params.Jsc = abs(currents(idx_sc));
            
            % 计算开路电压 (V)
            % 使用插值找到J=0的电压
            if any(currents > 0) && any(currents < 0)
                params.Voc = interp1(currents, voltages, 0);
            else
                % 如果电流没有变号，使用最大电压
                params.Voc = max(voltages);
            end
            
            % 计算功率密度 (mW/cm^2)
            power = voltages .* currents;
            
            % 找到最大功率点
            [params.Pmax, idx_mpp] = max(power);
            params.Vmpp = voltages(idx_mpp);
            params.Jmpp = currents(idx_mpp);
            
            % 计算填充因子
            if params.Jsc > 0 && params.Voc > 0
                params.FF = params.Pmax / (params.Jsc * params.Voc);
            else
                params.FF = 0;
            end
            
            % 计算光电转换效率 (%)
            params.PCE = 100 * params.Pmax / obj.AM15_power;
        end
        
        function jv_results = generateHysteresisJVCurve(obj, V_start, V_end, num_points, scan_rate)
            % 生成带迟滞的J-V特性曲线 (正向和反向扫描)
            % V_start: 起始电压 [V]
            % V_end: 终止电压 [V]
            % num_points: 电压点数
            % scan_rate: 扫描速率 [V/s]
            
            if nargin < 5
                scan_rate = 0.1; % 默认扫描速率
            end
            
            if nargin < 4
                num_points = 30; % 默认点数
            end
            
            if nargin < 3
                V_start = -0.2; % 默认起始电压
                V_end = 1.2;    % 默认终止电压
            end
            
            % 确保光照开启
            illumination_state = obj.params.illumination;
            obj.params.setIllumination(true);
            
            % 创建电压数组 (正向和反向扫描)
            voltages_forward = linspace(V_start, V_end, num_points);
            voltages_reverse = linspace(V_end, V_start, num_points);
            
            % 初始化电流密度数组
            currents_forward = zeros(size(voltages_forward));
            currents_reverse = zeros(size(voltages_reverse));
            
            % 计算每个电压点的电流密度 (正向扫描)
            disp('生成正向扫描J-V曲线...');
            for i = 1:length(voltages_forward)
                fprintf('计算电压点 %d/%d: %.3f V\n', i, length(voltages_forward), voltages_forward(i));
                
                % 设置施加电压
                obj.params.setAppliedVoltage(voltages_forward(i));
                
                % 求解漂移扩散方程
                results = obj.solver.solve();
                
                % 提取电流密度 (mA/cm^2)
                J_total = results.J_total(end, :);
                currents_forward(i) = J_total(1) * 1e3;
                
                % 模拟扫描速率
                if i < length(voltages_forward)
                    pause((voltages_forward(i+1) - voltages_forward(i)) / scan_rate);
                end
            end
            
            % 计算每个电压点的电流密度 (反向扫描)
            disp('生成反向扫描J-V曲线...');
            for i = 1:length(voltages_reverse)
                fprintf('计算电压点 %d/%d: %.3f V\n', i, length(voltages_reverse), voltages_reverse(i));
                
                % 设置施加电压
                obj.params.setAppliedVoltage(voltages_reverse(i));
                
                % 求解漂移扩散方程
                results = obj.solver.solve();
                
                % 提取电流密度 (mA/cm^2)
                J_total = results.J_total(end, :);
                currents_reverse(i) = J_total(1) * 1e3;
                
                % 模拟扫描速率
                if i < length(voltages_reverse)
                    pause((voltages_reverse(i) - voltages_reverse(i+1)) / scan_rate);
                end
            end
            
            % 恢复原始光照状态
            obj.params.setIllumination(illumination_state);
            
            % 计算正向扫描的性能参数
            jv_results_forward = obj.calculatePerformanceParameters(voltages_forward, currents_forward);
            
            % 计算反向扫描的性能参数
            jv_results_reverse = obj.calculatePerformanceParameters(voltages_reverse, currents_reverse);
            
            % 合并结果
            jv_results = struct();
            jv_results.forward = jv_results_forward;
            jv_results.reverse = jv_results_reverse;
            
            % 存储J-V数据
            jv_results.voltages_forward = voltages_forward;
            jv_results.currents_forward = currents_forward;
            jv_results.voltages_reverse = voltages_reverse;
            jv_results.currents_reverse = currents_reverse;
            
            % 计算迟滞指数 (HI = |PCE_reverse - PCE_forward| / PCE_reverse)
            jv_results.HI = abs(jv_results.reverse.PCE - jv_results.forward.PCE) / jv_results.reverse.PCE;
        end
    end
end
