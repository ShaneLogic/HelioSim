% test_script.m
% 简单测试脚本，用于验证主代码是否能正常运行

try
    % 运行主脚本
    disp('开始运行主脚本...');
    main_perovskite_cell;
    disp('主脚本运行成功！');
catch ME
    % 如果出错，显示错误信息
    disp('运行出错:');
    disp(ME.message);
    disp('错误位置:');
    disp(ME.stack(1));
end
