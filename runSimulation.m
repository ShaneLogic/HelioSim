%% runSimulation.m
% 太阳能电池模拟的主界面脚本，提供多种模拟选项

%% 判断是否使用图形界面
% 创建全局参数对象
params = SolarCellParams();

% 确定是否使用GUI界面
useGUI = true;
try
    % 检查是否可以使用GUIDE或App Designer
    v = ver('matlab');
    if str2double(v.Version) < 9.0
        useGUI = false;
    end
catch
    useGUI = false;
end

%% 根据界面类型运行模拟
if useGUI
    % 使用GUI界面
    solarCellGUI();
else
    % 使用命令行界面
    runCommandLine();
end

%% 命令行界面函数
function runCommandLine()
    % 创建全局参数对象
    params = SolarCellParams();
    
    % 显示菜单
    simType = menu('选择模拟类型:', ...
                   '平衡态', ...
                   '施加偏压', ...
                   '光照', ...
                   'J-V曲线', ...
                   '参数变化', ...
                   '退出');

    switch simType
        case 1 % 平衡态
            % 配置模拟
            config = struct('t_start', 0, ...
                            't_end', 1e-9, ...
                            'num_time_steps', 51, ...
                            'rel_tol', 1e-6, ...
                            'abs_tol', 1e-8, ...
                            'illumination', false, ...
                            'voltage_sweep', false);
            
            % 设置电压为平衡态(0V)
            params.setAppliedVoltage(0);
            
            % 创建求解器并运行模拟
            solver = DDSolver(params, config);
            results = solver.solve();
            
            % 可视化结果
            visualizer = Visualizer(params);
            visualizer.plotSteadyState(results);
            visualizer.plotBandDiagram(results);
            
        case 2 % 施加偏压
            % 从用户获取电压
            voltage = inputdlg('输入施加电压 (V):', '施加偏压', 1, {'0.5'});
            voltage = str2double(voltage{1});
            
            % 配置模拟
            config = struct('t_start', 0, ...
                            't_end', 1e-9, ...
                            'num_time_steps', 51, ...
                            'rel_tol', 1e-6, ...
                            'abs_tol', 1e-8, ...
                            'illumination', false, ...
                            'voltage_sweep', false);
            
            % 设置施加电压
            params.setAppliedVoltage(voltage);
            
            % 创建求解器并运行模拟
            solver = DDSolver(params, config);
            results = solver.solve();
            
            % 可视化结果
            visualizer = Visualizer(params);
            visualizer.plotSteadyState(results);
            visualizer.plotBandDiagram(results);
            visualizer.plotElectricField(results);
            visualizer.plotCurrentDensities(results);
            
        case 3 % 光照
            % 配置模拟
            config = struct('t_start', 0, ...
                            't_end', 1e-9, ...
                            'num_time_steps', 51, ...
                            'rel_tol', 1e-6, ...
                            'abs_tol', 1e-8, ...
                            'illumination', true, ...
                            'voltage_sweep', false);
            
            % 设置施加电压
            voltage = inputdlg('输入施加电压 (V):', '施加偏压', 1, {'0'});
            voltage = str2double(voltage{1});
            params.setAppliedVoltage(voltage);
            
            % 创建求解器
            solver = DDSolver(params, config);
            
            % 先运行暗态模拟
            params.setIllumination(false);
            dark_results = solver.solve();
            
            % 再运行光照模拟
            params.setIllumination(true);
            light_results = solver.solve();
            
            % 可视化结果
            visualizer = Visualizer(params);
            visualizer.plotComparisonWithLight(dark_results, light_results);
            
        case 4 % J-V曲线
            % 配置模拟
            config = struct('t_start', 0, ...
                            't_end', 1e-9, ...
                            'num_time_steps', 51, ...
                            'rel_tol', 1e-6, ...
                            'abs_tol', 1e-8, ...
                            'illumination', true, ...
                            'voltage_sweep', true);
            
            % 创建求解器
            solver = DDSolver(params, config);
            
            % 创建J-V分析器
            jv_analyzer = JVAnalyzer(params, solver);
            
            % 从用户获取电压范围
            prompt = {'起始电压 (V):', '结束电压 (V):', '数据点数:'};
            defaults = {'-0.2', '1.2', '20'};
            answer = inputdlg(prompt, 'J-V曲线设置', 1, defaults);
            
            if ~isempty(answer)
                V_start = str2double(answer{1});
                V_end = str2double(answer{2});
                V_points = round(str2double(answer{3}));
                
                % 生成J-V曲线
                jv_results = jv_analyzer.generateJVCurve(V_start, V_end, V_points);
            end
            
        case 5 % 参数变化
            % 选择要变化的参数
            paramType = menu('选择要变化的参数:', ...
                             '吸收层厚度', ...
                             '吸收层带隙', ...
                             'ETL掺杂', ...
                             'HTL掺杂', ...
                             '迁移率', ...
                             '返回主菜单');
            
            if paramType == 6
                runSimulation;
                return;
            end
            
            % 配置模拟
            config = struct('t_start', 0, ...
                            't_end', 1e-9, ...
                            'num_time_steps', 51, ...
                            'rel_tol', 1e-6, ...
                            'abs_tol', 1e-8, ...
                            'illumination', true, ...
                            'voltage_sweep', false);
            
            % 获取变化次数
            n_var = inputdlg('输入变化次数:', '参数变化', 1, {'3'});
            n_var = round(str2double(n_var{1}));
            
            % 初始化结果数组
            Voc_array = zeros(1, n_var);
            Jsc_array = zeros(1, n_var);
            FF_array = zeros(1, n_var);
            PCE_array = zeros(1, n_var);
            param_values = zeros(1, n_var);
            
            % 创建求解器
            solver = DDSolver(params, config);
            
            % 创建J-V分析器
            jv_analyzer = JVAnalyzer(params, solver);
            
            % 循环变化参数
            for i = 1:n_var
                % 根据参数类型设置参数值
                switch paramType
                    case 1 % 吸收层厚度
                        % 从200 nm到700 nm
                        thickness = 2e-5 + (i-1) * 5e-5/(n_var-1);
                        params.setAbsorberThickness(thickness);
                        param_values(i) = thickness * 1e4; % 转换为µm
                        param_name = '吸收层厚度 (µm)';
                        
                    case 2 % 吸收层带隙
                        % 从1.2 eV到2.0 eV
                        bandgap = 1.2 + (i-1) * 0.8/(n_var-1);
                        params.setAbsorberBandgap(bandgap);
                        param_values(i) = bandgap;
                        param_name = '吸收层带隙 (eV)';
                        
                    case 3 % ETL掺杂
                        % 从1e16到1e18 cm^-3
                        doping = 1e16 * 10^((i-1) * 2/(n_var-1));
                        params.setETLDoping(doping);
                        param_values(i) = log10(doping);
                        param_name = 'ETL掺杂浓度 log_{10}(cm^{-3})';
                        
                    case 4 % HTL掺杂
                        % 从1e16到1e18 cm^-3
                        doping = 1e16 * 10^((i-1) * 2/(n_var-1));
                        params.setHTLDoping(doping);
                        param_values(i) = log10(doping);
                        param_name = 'HTL掺杂浓度 log_{10}(cm^{-3})';
                        
                    case 5 % 迁移率
                        % 从0.1到10 cm^2/Vs
                        mobility = 0.1 * 10^((i-1) * 2/(n_var-1));
                        params.setAbsorberMobility(mobility);
                        param_values(i) = mobility;
                        param_name = '吸收层迁移率 (cm^2/Vs)';
                end
                
                % 设置求解器
                solver = DDSolver(params, config);
                jv_analyzer = JVAnalyzer(params, solver);
                
                % 生成J-V曲线
                jv_results = jv_analyzer.generateJVCurve(-0.2, 1.2, 10);
                
                % 存储结果
                Voc_array(i) = jv_results.Voc;
                Jsc_array(i) = abs(jv_results.Jsc);
                FF_array(i) = jv_results.FF;
                PCE_array(i) = jv_results.PCE;
            end
            
            % 绘制参数变化结果
            figure('Name', '参数变化结果', 'Position', [100, 100, 1200, 800]);
            
            % 绘制Voc
            subplot(2,2,1);
            plot(param_values, Voc_array, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
            grid on;
            xlabel(param_name, 'FontSize', 12);
            ylabel('V_{oc} (V)', 'FontSize', 12);
            title('开路电压', 'FontSize', 14);
            
            % 绘制Jsc
            subplot(2,2,2);
            plot(param_values, Jsc_array, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
            grid on;
            xlabel(param_name, 'FontSize', 12);
            ylabel('J_{sc} (mA/cm^2)', 'FontSize', 12);
            title('短路电流密度', 'FontSize', 14);
            
            % 绘制FF
            subplot(2,2,3);
            plot(param_values, FF_array*100, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
            grid on;
            xlabel(param_name, 'FontSize', 12);
            ylabel('填充因子 (%)', 'FontSize', 12);
            title('填充因子', 'FontSize', 14);
            
            % 绘制PCE
            subplot(2,2,4);
            plot(param_values, PCE_array, 'ko-', 'LineWidth', 2, 'MarkerSize', 8);
            grid on;
            xlabel(param_name, 'FontSize', 12);
            ylabel('PCE (%)', 'FontSize', 12);
            title('光电转换效率', 'FontSize', 14);
            
        case 6 % 退出
            close all;
            clc;
            return;
    end
end

%% GUI界面函数
function solarCellGUI()
    % 创建全局参数对象
    params = SolarCellParams();
    
    % 创建图形界面的主窗口
    mainFig = figure('Name', 'HelioSim - 太阳能电池模拟', ...
                    'NumberTitle', 'off', ...
                    'Position', [100, 100, 800, 600], ...
                    'MenuBar', 'none', ...
                    'ToolBar', 'none');
                
    % 创建模拟类型选择面板
    simPanel = uipanel('Title', '模拟类型', ...
                    'Position', [0.05, 0.75, 0.9, 0.2]);
                
    % 创建模拟类型选择按钮
    btnEquil = uicontrol('Parent', simPanel, ...
                        'Style', 'pushbutton', ...
                        'String', '平衡态模拟', ...
                        'Position', [30, 50, 120, 40], ...
                        'Callback', @runEquilibriumSim);
                    
    btnBias = uicontrol('Parent', simPanel, ...
                        'Style', 'pushbutton', ...
                        'String', '施加偏压', ...
                        'Position', [180, 50, 120, 40], ...
                        'Callback', @runBiasSim);
                    
    btnLight = uicontrol('Parent', simPanel, ...
                        'Style', 'pushbutton', ...
                        'String', '光照模拟', ...
                        'Position', [330, 50, 120, 40], ...
                        'Callback', @runIlluminationSim);
                    
    btnJV = uicontrol('Parent', simPanel, ...
                    'Style', 'pushbutton', ...
                    'String', 'J-V曲线', ...
                    'Position', [480, 50, 120, 40], ...
                    'Callback', @runJVSim);
                
    btnParam = uicontrol('Parent', simPanel, ...
                        'Style', 'pushbutton', ...
                        'String', '参数变化', ...
                        'Position', [630, 50, 120, 40], ...
                        'Callback', @runParamVarSim);
                    
    % 创建材料参数设置面板
    paramPanel = uipanel('Title', '材料参数设置', ...
                        'Position', [0.05, 0.35, 0.9, 0.35]);
                    
    % 吸收层参数
    absText = uicontrol('Parent', paramPanel, ...
                        'Style', 'text', ...
                        'String', '吸收层参数', ...
                        'Position', [30, 120, 100, 20], ...
                        'HorizontalAlignment', 'left');
                    
    uicontrol('Parent', paramPanel, ...
            'Style', 'text', ...
            'String', '厚度 (nm):', ...
            'Position', [30, 90, 80, 20], ...
            'HorizontalAlignment', 'left');
        
    thickEdit = uicontrol('Parent', paramPanel, ...
                        'Style', 'edit', ...
                        'String', '500', ...
                        'Position', [120, 90, 60, 20]);
                    
    uicontrol('Parent', paramPanel, ...
            'Style', 'text', ...
            'String', '带隙 (eV):', ...
            'Position', [30, 60, 80, 20], ...
            'HorizontalAlignment', 'left');
        
    bandgapEdit = uicontrol('Parent', paramPanel, ...
                            'Style', 'edit', ...
                            'String', '1.5', ...
                            'Position', [120, 60, 60, 20]);
                        
    uicontrol('Parent', paramPanel, ...
            'Style', 'text', ...
            'String', '迁移率 (cm²/Vs):', ...
            'Position', [30, 30, 100, 20], ...
            'HorizontalAlignment', 'left');
        
    mobilityEdit = uicontrol('Parent', paramPanel, ...
                            'Style', 'edit', ...
                            'String', '10', ...
                            'Position', [140, 30, 60, 20]);
                        
    % ETL参数
    etlText = uicontrol('Parent', paramPanel, ...
                        'Style', 'text', ...
                        'String', 'ETL参数', ...
                        'Position', [250, 120, 100, 20], ...
                        'HorizontalAlignment', 'left');
                    
    uicontrol('Parent', paramPanel, ...
            'Style', 'text', ...
            'String', '掺杂浓度 (cm⁻³):', ...
            'Position', [250, 90, 120, 20], ...
            'HorizontalAlignment', 'left');
        
    etlDopingEdit = uicontrol('Parent', paramPanel, ...
                            'Style', 'edit', ...
                            'String', '1e17', ...
                            'Position', [380, 90, 60, 20]);
                        
    % HTL参数
    htlText = uicontrol('Parent', paramPanel, ...
                        'Style', 'text', ...
                        'String', 'HTL参数', ...
                        'Position', [470, 120, 100, 20], ...
                        'HorizontalAlignment', 'left');
                    
    uicontrol('Parent', paramPanel, ...
            'Style', 'text', ...
            'String', '掺杂浓度 (cm⁻³):', ...
            'Position', [470, 90, 120, 20], ...
            'HorizontalAlignment', 'left');
        
    htlDopingEdit = uicontrol('Parent', paramPanel, ...
                            'Style', 'edit', ...
                            'String', '1e17', ...
                            'Position', [600, 90, 60, 20]);
                        
    % 应用参数按钮
    applyBtn = uicontrol('Style', 'pushbutton', ...
                        'String', '应用参数', ...
                        'Position', [350, 150, 100, 30], ...
                        'Callback', @applyParameters);
                    
    % 状态显示区域
    statusPanel = uipanel('Title', '状态信息', ...
                        'Position', [0.05, 0.05, 0.9, 0.25]);
                    
    statusText = uicontrol('Parent', statusPanel, ...
                            'Style', 'text', ...
                            'String', '准备就绪，请选择模拟类型', ...
                            'Position', [20, 50, 700, 30], ...
                            'HorizontalAlignment', 'left');
    
    % 模拟类型回调函数
    function runEquilibriumSim(~, ~)
        % 更新状态
        set(statusText, 'String', '正在运行平衡态模拟...');
        drawnow;
        
        % 配置模拟
        config = struct('t_start', 0, ...
                       't_end', 1e-9, ...
                       'num_time_steps', 51, ...
                       'rel_tol', 1e-6, ...
                       'abs_tol', 1e-8, ...
                       'illumination', false, ...
                       'voltage_sweep', false);
        
        % 设置电压为0V（平衡态）
        params.setAppliedVoltage(0);
        
        % 创建求解器并运行模拟
        solver = DDSolver(params, config);
        try
            results = solver.solve();
            
            % 可视化结果
            visualizer = Visualizer(params);
            visualizer.plotSteadyState(results);
            visualizer.plotBandDiagram(results);
            
            % 更新状态
            set(statusText, 'String', '平衡态模拟完成');
        catch ME
            set(statusText, 'String', ['模拟错误: ' ME.message]);
        end
    end

    function runBiasSim(~, ~)
        % 创建电压输入对话框
        voltageInput = inputdlg('输入施加电压 (V):', '施加偏压', 1, {'0.5'});
        if isempty(voltageInput)
            return;
        end
        
        voltage = str2double(voltageInput{1});
        
        % 更新状态
        set(statusText, 'String', ['正在运行施加偏压 ' num2str(voltage) 'V 的模拟...']);
        drawnow;
        
        % 配置模拟
        config = struct('t_start', 0, ...
                       't_end', 1e-9, ...
                       'num_time_steps', 51, ...
                       'rel_tol', 1e-6, ...
                       'abs_tol', 1e-8, ...
                       'illumination', false, ...
                       'voltage_sweep', false);
        
        % 设置施加电压
        params.setAppliedVoltage(voltage);
        
        % 创建求解器并运行模拟
        solver = DDSolver(params, config);
        try
            results = solver.solve();
            
            % 可视化结果
            visualizer = Visualizer(params);
            visualizer.plotSteadyState(results);
            visualizer.plotBandDiagram(results);
            visualizer.plotElectricField(results);
            visualizer.plotCurrentDensities(results);
            
            % 更新状态
            set(statusText, 'String', ['施加偏压 ' num2str(voltage) 'V 的模拟完成']);
        catch ME
            set(statusText, 'String', ['模拟错误: ' ME.message]);
        end
    end

    function runIlluminationSim(~, ~)
        % 创建电压输入对话框
        voltageInput = inputdlg('输入施加电压 (V):', '光照模拟', 1, {'0'});
        if isempty(voltageInput)
            return;
        end
        
        voltage = str2double(voltageInput{1});
        
        % 更新状态
        set(statusText, 'String', ['正在运行光照模拟(电压=' num2str(voltage) 'V)...']);
        drawnow;
        
        % 配置模拟
        config = struct('t_start', 0, ...
                       't_end', 1e-9, ...
                       'num_time_steps', 51, ...
                       'rel_tol', 1e-6, ...
                       'abs_tol', 1e-8, ...
                       'illumination', true, ...
                       'voltage_sweep', false);
        
        % 设置施加电压
        params.setAppliedVoltage(voltage);
        
        % 创建求解器
        solver = DDSolver(params, config);
        
        try
            % 先运行暗态模拟
            params.setIllumination(false);
            dark_results = solver.solve();
            
            % 再运行光照模拟
            params.setIllumination(true);
            light_results = solver.solve();
            
            % 可视化结果
            visualizer = Visualizer(params);
            visualizer.plotComparisonWithLight(dark_results, light_results);
            
            % 更新状态
            set(statusText, 'String', '光照模拟完成');
        catch ME
            set(statusText, 'String', ['模拟错误: ' ME.message]);
        end
    end

    function runJVSim(~, ~)
        % 创建J-V设置对话框
        prompt = {'起始电压 (V):', '结束电压 (V):', '数据点数:'};
        defaults = {'-0.2', '1.2', '20'};
        jvSettings = inputdlg(prompt, 'J-V曲线设置', 1, defaults);
        
        if isempty(jvSettings)
            return;
        end
        
        V_start = str2double(jvSettings{1});
        V_end = str2double(jvSettings{2});
        V_points = round(str2double(jvSettings{3}));
        
        % 更新状态
        set(statusText, 'String', '正在生成J-V曲线...');
        drawnow;
        
        % 配置模拟
        config = struct('t_start', 0, ...
                       't_end', 1e-9, ...
                       'num_time_steps', 51, ...
                       'rel_tol', 1e-6, ...
                       'abs_tol', 1e-8, ...
                       'illumination', true, ...
                       'voltage_sweep', true);
        
        % 创建求解器
        solver = DDSolver(params, config);
        
        % 创建J-V分析器
        jv_analyzer = JVAnalyzer(params, solver);
        
        try
            % 生成J-V曲线
            jv_results = jv_analyzer.generateJVCurve(V_start, V_end, V_points);
            
            % 更新状态
            set(statusText, 'String', ['J-V曲线生成完成，Voc = ' num2str(jv_results.Voc, '%.3f') ' V, ' ...
                                     'Jsc = ' num2str(abs(jv_results.Jsc), '%.2f') ' mA/cm², ' ...
                                     'FF = ' num2str(jv_results.FF*100, '%.1f') ' %, ' ...
                                     'PCE = ' num2str(jv_results.PCE, '%.2f') ' %']);
        catch ME
            set(statusText, 'String', ['模拟错误: ' ME.message]);
        end
    end

    function runParamVarSim(~, ~)
        % 创建参数变化类型选择对话框
        paramTypes = {'吸收层厚度', '吸收层带隙', 'ETL掺杂', 'HTL掺杂', '迁移率'};
        [paramTypeIdx, ok] = listdlg('ListString', paramTypes, ...
                                    'SelectionMode', 'single', ...
                                    'Name', '选择参数变化类型', ...
                                    'PromptString', '选择要变化的参数:');
        
        if ~ok || isempty(paramTypeIdx)
            return;
        end
        
        % 获取变化次数
        varCountInput = inputdlg('输入变化次数:', '参数变化', 1, {'3'});
        if isempty(varCountInput)
            return;
        end
        
        n_var = round(str2double(varCountInput{1}));
        
        % 更新状态
        set(statusText, 'String', ['正在进行' paramTypes{paramTypeIdx} '参数变化模拟...']);
        drawnow;
        
        % 配置模拟
        config = struct('t_start', 0, ...
                       't_end', 1e-9, ...
                       'num_time_steps', 51, ...
                       'rel_tol', 1e-6, ...
                       'abs_tol', 1e-8, ...
                       'illumination', true, ...
                       'voltage_sweep', false);
        
        % 初始化结果数组
        Voc_array = zeros(1, n_var);
        Jsc_array = zeros(1, n_var);
        FF_array = zeros(1, n_var);
        PCE_array = zeros(1, n_var);
        param_values = zeros(1, n_var);
        
        % 创建求解器
        solver = DDSolver(params, config);
        
        % 创建J-V分析器
        jv_analyzer = JVAnalyzer(params, solver);
        
        try
            % 循环变化参数
            for i = 1:n_var
                % 根据参数类型设置参数值
                switch paramTypeIdx
                    case 1 % 吸收层厚度
                        % 从200 nm到700 nm
                        thickness = 2e-5 + (i-1) * 5e-5/(n_var-1);
                        params.setAbsorberThickness(thickness);
                        param_values(i) = thickness * 1e4; % 转换为µm
                        param_name = '吸收层厚度 (µm)';
                        
                    case 2 % 吸收层带隙
                        % 从1.2 eV到2.0 eV
                        bandgap = 1.2 + (i-1) * 0.8/(n_var-1);
                        params.setAbsorberBandgap(bandgap);
                        param_values(i) = bandgap;
                        param_name = '吸收层带隙 (eV)';
                        
                    case 3 % ETL掺杂
                        % 从1e16到1e18 cm^-3
                        doping = 1e16 * 10^((i-1) * 2/(n_var-1));
                        params.setETLDoping(doping);
                        param_values(i) = log10(doping);
                        param_name = 'ETL掺杂浓度 log_{10}(cm^{-3})';
                        
                    case 4 % HTL掺杂
                        % 从1e16到1e18 cm^-3
                        doping = 1e16 * 10^((i-1) * 2/(n_var-1));
                        params.setHTLDoping(doping);
                        param_values(i) = log10(doping);
                        param_name = 'HTL掺杂浓度 log_{10}(cm^{-3})';
                        
                    case 5 % 迁移率
                        % 从0.1到10 cm^2/Vs
                        mobility = 0.1 * 10^((i-1) * 2/(n_var-1));
                        params.setAbsorberMobility(mobility);
                        param_values(i) = mobility;
                        param_name = '吸收层迁移率 (cm^2/Vs)';
                end
                
                % 更新状态
                set(statusText, 'String', ['参数变化模拟进度: ' num2str(i) '/' num2str(n_var)]);
                drawnow;
                
                % 生成J-V曲线
                jv_results = jv_analyzer.generateJVCurve(-0.2, 1.2, 10);
                
                % 存储结果
                Voc_array(i) = jv_results.Voc;
                Jsc_array(i) = abs(jv_results.Jsc);
                FF_array(i) = jv_results.FF;
                PCE_array(i) = jv_results.PCE;
            end
            
            % 绘制参数变化结果
            paramVarFig = figure('Name', '参数变化结果', 'Position', [100, 100, 1200, 800]);
            
            % 绘制Voc
            subplot(2,2,1);
            plot(param_values, Voc_array, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
            grid on;
            xlabel(param_name, 'FontSize', 12);
            ylabel('V_{oc} (V)', 'FontSize', 12);
            title('开路电压', 'FontSize', 14);
            
            % 绘制Jsc
            subplot(2,2,2);
            plot(param_values, Jsc_array, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
            grid on;
            xlabel(param_name, 'FontSize', 12);
            ylabel('J_{sc} (mA/cm^2)', 'FontSize', 12);
            title('短路电流密度', 'FontSize', 14);
            
            % 绘制FF
            subplot(2,2,3);
            plot(param_values, FF_array*100, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
            grid on;
            xlabel(param_name, 'FontSize', 12);
            ylabel('填充因子 (%)', 'FontSize', 12);
            title('填充因子', 'FontSize', 14);
            
            % 绘制PCE
            subplot(2,2,4);
            plot(param_values, PCE_array, 'ko-', 'LineWidth', 2, 'MarkerSize', 8);
            grid on;
            xlabel(param_name, 'FontSize', 12);
            ylabel('PCE (%)', 'FontSize', 12);
            title('光电转换效率', 'FontSize', 14);
            
            % 更新状态
            set(statusText, 'String', ['参数变化模拟完成: ' param_name]);
        catch ME
            set(statusText, 'String', ['模拟错误: ' ME.message]);
        end
    end

    function applyParameters(~, ~)
        % 获取GUI中的参数值
        try
            % 吸收层参数
            thickness = str2double(get(thickEdit, 'String')) * 1e-7; % 从nm转为cm
            bandgap = str2double(get(bandgapEdit, 'String'));
            mobility = str2double(get(mobilityEdit, 'String'));
            
            % ETL参数
            etlDoping = str2double(get(etlDopingEdit, 'String'));
            
            % HTL参数
            htlDoping = str2double(get(htlDopingEdit, 'String'));
            
            % 设置参数
            params.setAbsorberThickness(thickness);
            params.setAbsorberBandgap(bandgap);
            params.setAbsorberMobility(mobility);
            params.setETLDoping(etlDoping);
            params.setHTLDoping(htlDoping);
            
            % 更新状态
            set(statusText, 'String', '参数已应用');
        catch ME
            set(statusText, 'String', ['参数设置错误: ' ME.message]);
        end
    end
end