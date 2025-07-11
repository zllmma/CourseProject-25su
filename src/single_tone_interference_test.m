% =========================================================================
%               单频干扰测试：七种频率估计算法的抗干扰性能分析
%               Single-Tone Interference Test for Frequency Estimation
% =========================================================================
%
% 功能:
%   分析七种频率估计算法在不同干扰频率和不同窗函数下的性能表现
%
% 使用示例:
%   window_type = 'hann';              % 汉宁窗
%   ISR_dB = 0;                        % 干信比 (干扰功率/信号功率)
% =========================================================================

clear;
close all;
clc;

addpath('algorithms');  % 添加算法目录到 MATLAB 路径

% --- 0. 用户参数选择 (修改此处) ---
window_type = 'hann';              % 可选：'none', 'hann', 'hamming', 'blackman', 'kaiser'
ISR_dB = 0;                        % 干信比 (干扰功率/信号功率, dB)

% --- 1. 仿真参数设置 ---
fs = 200e6;              % 采样频率 (Hz)
N = 1024;                % 采样点数
t = (0:N-1) / fs;        % 时间向量

% 设置一个非FFT整数倍的频率，以突显栅栏效应
f_center = 50e6; % 中心频率 (50 MHz)
delta_f0 = fs / N; % 频率分辨率 (Hz)
offset = 0.1; % 相对频偏
f_true = f_center + offset * delta_f0; % 真实信号频率 (Hz)

f_start = 49e6;          % 干扰起始频率 (Hz)
f_end = 51e6;            % 干扰结束频率 (Hz)
f_step = 0.2e6;           % 干扰频率步长 (200kHz)
f_interference_list = f_start:f_step:f_end; % 干扰频率列表

num_interferences = length(f_interference_list);
num_trials = 1000;        % 蒙特卡洛试验次数 (减少计算量)

% CZT 和 改进CZT 算法的参数
q = 1;                   % 频率细化区间大小控制参数
M = 64;                  % 频率细化倍数

% --- 2. 生成窗函数并进行归一化 ---
switch window_type
    case 'none'
        w = ones(1, N);  % 不使用窗函数，相当于矩形窗但不归一化
        w_name = '无窗函数';
    case 'hann'
        w = hann(N)';
        w_name = 'Hann窗';
    case 'hamming'
        w = hamming(N)';
        w_name = 'Hamming窗';
    case 'blackman'
        w = blackman(N)';
        w_name = 'Blackman窗';
    case 'kaiser'
        beta = 8.6;  % Kaiser窗参数
        w = kaiser(N, beta)';
        w_name = 'Kaiser窗';
    otherwise
        error('未知窗类型: %s。支持的窗类型: none, hann, hamming, blackman, kaiser', window_type);
end

% 窗函数归一化 (保持信号功率不变)
% 对于'none'选项不进行归一化，保持原始信号的幅度
if ~strcmp(window_type, 'none')
    w = w / sqrt(mean(w.^2));
end

% 计算干信比对应的幅度因子
A_interference = 10^(ISR_dB/20);  % 干扰信号幅度 (信号幅度=1)

% 初始化结果存储向量 - 每个算法一个向量
rmse_fft = zeros(1, num_interferences);
rmse_czt = zeros(1, num_interferences);
rmse_improved_czt = zeros(1, num_interferences);
rmse_rife = zeros(1, num_interferences);
rmse_mrife = zeros(1, num_interferences);
rmse_irife = zeros(1, num_interferences);
rmse_iirife = zeros(1, num_interferences);

% --- 3. 执行蒙特卡洛模拟 ---
fprintf('开始蒙特卡洛模拟，窗类型: %s，干信比: %d dB...\n', w_name, ISR_dB);

parfor idx = 1:num_interferences
    fprintf('处理干扰频率: %.2f MHz\n', f_interference_list(idx)/1e6);
    f_interference = f_interference_list(idx);
    
    % 用于存储当前干扰频率下各算法的误差
    errors_fft = zeros(1, num_trials);
    errors_czt = zeros(1, num_trials);
    errors_improved_czt = zeros(1, num_trials);
    errors_rife = zeros(1, num_trials);
    errors_mrife = zeros(1, num_trials);
    errors_irife = zeros(1, num_trials);
    errors_iirife = zeros(1, num_trials);
    
    phases = 2 * pi * rand(1, num_trials);        % 信号随机相位
    interf_phases = 2 * pi * rand(1, num_trials); % 干扰随机相位
    
    for j = 1:num_trials
        % a. 生成纯净信号 (使用复正弦信号)
        phi = phases(j);
        s_clean = exp(1j * (2 * pi * f_true * t + phi));
        
        % b. 生成单频干扰 (固定幅度)
        phi_interf = interf_phases(j);
        interference = A_interference * exp(1j * (2 * pi * f_interference * t + phi_interf));
        
        % c. 合成信号 (无背景噪声)
        s_noisy = s_clean + interference;
        
        % d. 应用窗函数
        s_windowed = s_noisy .* w;
        
        % e. 使用七种算法进行频率估计
        f_fft = fft_est(s_windowed, fs);
        f_czt = czt_est(s_windowed, fs, q, M);
        f_improved_czt = improved_czt_est(s_windowed, fs, q, M);
        f_rife = rife_est(s_windowed, fs);
        f_mrife = mrife_est(s_windowed, fs);
        f_irife = irife_est(s_windowed, fs);
        f_iirife = iirife_est(s_windowed, fs);
        
        % f. 计算并存储误差
        errors_fft(j) = f_fft - f_true;
        errors_czt(j) = f_czt - f_true;
        errors_improved_czt(j) = f_improved_czt - f_true;
        errors_rife(j) = f_rife - f_true;
        errors_mrife(j) = f_mrife - f_true;
        errors_irife(j) = f_irife - f_true;
        errors_iirife(j) = f_iirife - f_true;
    end
    
    % g. 计算当前干扰频率下的RMSE
    rmse_fft(idx) = sqrt(mean(errors_fft.^2));
    rmse_czt(idx) = sqrt(mean(errors_czt.^2));
    rmse_improved_czt(idx) = sqrt(mean(errors_improved_czt.^2));
    rmse_rife(idx) = sqrt(mean(errors_rife.^2));
    rmse_mrife(idx) = sqrt(mean(errors_mrife.^2));
    rmse_irife(idx) = sqrt(mean(errors_irife.^2));
    rmse_iirife(idx) = sqrt(mean(errors_iirife.^2));
end

fprintf('模拟完成。\n');

% --- 4. 绘制结果 ---
% 在一个窗口中绘制RMSE结果
figure('Position', [100, 100, 1200, 700]);

% RMSE对比图
semilogy(f_interference_list/1e6, rmse_fft, '-o', 'LineWidth', 1.5, 'DisplayName', 'FFT-Peak');
hold on;
semilogy(f_interference_list/1e6, rmse_czt, '-s', 'LineWidth', 1.5, 'DisplayName', 'CZT');
semilogy(f_interference_list/1e6, rmse_improved_czt, '-^', 'LineWidth', 1.5, 'DisplayName', '改进 CZT');
semilogy(f_interference_list/1e6, rmse_rife, '-d', 'LineWidth', 1.5, 'DisplayName', 'RIFE');
semilogy(f_interference_list/1e6, rmse_mrife, '-x', 'LineWidth', 1.5, 'DisplayName', 'MRIFE');
semilogy(f_interference_list/1e6, rmse_irife, '-+', 'LineWidth', 1.5, 'DisplayName', 'IRIFE');
semilogy(f_interference_list/1e6, rmse_iirife, '-*', 'LineWidth', 1.5, 'DisplayName', 'IIRIFE');

% 添加垂直线标记信号频率
line([f_true/1e6, f_true/1e6], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% 添加信号频率标注
text(f_true/1e6 + 0.05, max(ylim)*0.8, sprintf('信号频率: %.4f MHz', f_true/1e6), ...
    'FontSize', 10, ...
    'BackgroundColor', 'white');

hold off;
grid on;
title(sprintf('频率估计算法在单频干扰下的性能 (频偏: %s, 窗函数: %s, 干信比: %d dB)', ...
    num2str(offset), w_name, ISR_dB));
xlabel('干扰频率 (MHz)');
ylabel('RMSE (Hz)');
legend('show', 'Location', 'best');

