% =========================================================================
%               七种频率估计算法的性能对比脚本 (支持多种干扰和窗函数)
% =========================================================================
%
% 新增功能:
%   1. 可选择干扰类型：'awgn'(高斯白噪声), 'sinusoidal'(单频干扰), 
%      'chirp'(线性调频干扰), 'impulse'(脉冲干扰)
%   2. 可选择窗函数：'none'(不使用窗), 'rect'(矩形窗), 'hann'(汉宁窗), 'hamming'(海明窗), 'blackman'(布莱克曼窗)
%
% 使用示例:
%   interference_type = 'sinusoidal';  % 单频干扰
%   window_type = 'hann';              % 汉宁窗
%   p_impulse = 0.1;                   % 脉冲干扰出现概率
%   f_interference = 55e6;             % 单频干扰频率 (Hz)
%   f0_chirp = 40e6;                   % 线性调频起始频率 (Hz)
%   k_chirp = 20e12;                   % 线性调频斜率 (Hz/s)
% =========================================================================

clear;
close all;
clc;

% 添加 algorithms 目录到 MATLAB 路径
addpath('algorithms');

% --- 0. 用户参数选择 (修改此处) ---
interference_type = 'chirp';  % 可选：'awgn', 'sinusoidal', 'chirp', 'impulse'
window_type = 'rect';              % 可选：'none', 'rect', 'hann', 'hamming', 'blackman'
p_impulse = 0.1;                   % 脉冲干扰出现概率 (仅interference_type='impulse'时有效)
f_interference = 55e6;             % 单频干扰频率 (Hz) (仅interference_type='sinusoidal'时有效)
f0_chirp = 40e6;                   % 线性调频起始频率 (Hz) (仅interference_type='chirp'时有效)
k_chirp = 20e12;                   % 线性调频斜率 (Hz/s) (仅interference_type='chirp'时有效)

% --- 1. 仿真参数设置 ---
fs = 200e6;              % 采样频率 (Hz)
N = 1024;                % 采样点数
t = (0:N-1) / fs;        % 时间向量

% 设置一个非FFT整数倍的频率，以突显栅栏效应
f_center = 50e6; % 中心频率 (50 MHz)
delta_f0 = fs / N; % 频率分辨率 (Hz)
offset = 0.4; % 相对频偏
f_true = f_center + offset * delta_f0; % 真实信号频率 (Hz)

SIR_dB = -20:2:20;       % 信干比范围 (dB)
num_trials = 1000;       % 每个SIR下的蒙特卡洛试验次数

% CZT 和 改进CZT 算法的参数
q = 1;                   % 频率细化区间大小控制参数
M = 64;                  % 频率细化倍数

% --- 2. 生成窗函数并进行归一化 ---
switch window_type
    case 'none'
        w = ones(1, N);  % 不使用窗函数，相当于矩形窗但不归一化
        w_name = '无';
    case 'rect'
        w = ones(1, N);
        w_name = '矩形窗';
    case 'hann'
        w = hann(N)';
        w_name = '汉宁窗';
    case 'hamming'
        w = hamming(N)';
        w_name = '海明窗';
    case 'blackman'
        w = blackman(N)';
        w_name = '布莱克曼窗';
    otherwise
        error('未知窗类型: %s', window_type);
end

% 窗函数归一化 (保持信号功率不变)
% 对于'none'选项不进行归一化，保持原始信号的幅度
if ~strcmp(window_type, 'none')
    w = w / sqrt(mean(w.^2));
end

% --- 3. 初始化结果存储变量 ---
num_sirs = length(SIR_dB);
rmse_fft = zeros(1, num_sirs);
rmse_czt = zeros(1, num_sirs);
rmse_improved_czt = zeros(1, num_sirs);
rmse_rife = zeros(1, num_sirs);
rmse_mrife = zeros(1, num_sirs);
rmse_irife = zeros(1, num_sirs);
rmse_iirife = zeros(1, num_sirs);

% 添加标准差存储变量
std_fft = zeros(1, num_sirs);
std_czt = zeros(1, num_sirs);
std_improved_czt = zeros(1, num_sirs);
std_rife = zeros(1, num_sirs);
std_mrife = zeros(1, num_sirs);
std_irife = zeros(1, num_sirs);
std_iirife = zeros(1, num_sirs);

% --- 4. 执行蒙特卡洛模拟 ---
fprintf('开始蒙特卡洛模拟，干扰类型: %s, 窗类型: %s...\n', interference_type, w_name);

% 并行执行每个SIR值的试验
parfor i = 1:num_sirs
    sir_current_db = SIR_dB(i);
    fprintf('正在处理 SIR = %.0f dB...\n', sir_current_db);
    
    % 用于存储单次SIR下的所有试验误差
    errors_fft = zeros(1, num_trials);
    errors_czt = zeros(1, num_trials);
    errors_improved_czt = zeros(1, num_trials);
    errors_rife = zeros(1, num_trials);
    errors_mrife = zeros(1, num_trials);
    errors_irife = zeros(1, num_trials);
    errors_iirife = zeros(1, num_trials);
    phases = 2 * pi * rand(1, num_trials); % 随机相位用于每次试验

    for j = 1:num_trials
        % a. 生成纯净信号 (使用复正弦信号)
        phi = phases(j);
        s_clean = exp(1j * (2 * pi * f_true * t + phi));
        
        % b. 根据SIR计算干扰功率并生成干扰
        signal_power = mean(abs(s_clean).^2);
        sir_linear = 10^(sir_current_db / 10);
        interference_power = signal_power / sir_linear;
        interference = zeros(1, N); % 初始化干扰信号
        switch interference_type
            case 'awgn'
                % 高斯白噪声
                interference_real = randn(1, N) * sqrt(interference_power/2);
                interference_imag = randn(1, N) * sqrt(interference_power/2);
                interference = interference_real + 1j * interference_imag;
                
            case 'sinusoidal'
                % 单频干扰 (带随机相位)
                phi_interf = 2*pi*rand();
                A = sqrt(2*interference_power); % 幅度计算
                interference = A * exp(1j * (2 * pi * f_interference * t + phi_interf));
                
            case 'chirp'
                % 线性调频干扰 (带随机相位)
                phi_interf = 2*pi*rand();
                A = sqrt(2*interference_power); % 幅度计算
                interference = A * exp(1j * (2 * pi * (f0_chirp * t + 0.5 * k_chirp * t.^2) + phi_interf));
                
            case 'impulse'
                % 脉冲干扰
                sigma_impulse = sqrt(interference_power/(2*p_impulse));
                impulse_real = double(rand(1, N) < p_impulse) .* ...
                              (sigma_impulse * randn(1, N));
                impulse_imag = double(rand(1, N) < p_impulse) .* ...
                              (sigma_impulse * randn(1, N));
                interference = impulse_real + 1j * impulse_imag;
        end
        
        s_noisy = s_clean + interference;
        
        % c. 应用窗函数
        s_windowed = s_noisy .* w;
        
        % d. 使用七种算法进行频率估计
        f_fft = fft_est(s_windowed, fs);
        f_czt = czt_est(s_windowed, fs, q, M);
        f_improved_czt = improved_czt_est(s_windowed, fs, q, M);
        f_rife = rife_est(s_windowed, fs);
        f_mrife = mrife_est(s_windowed, fs);
        f_irife = irife_est(s_windowed, fs);
        f_iirife = iirife_est(s_windowed, fs);

        % e. 计算并存储误差
        errors_fft(j) = f_fft - f_true;
        errors_czt(j) = f_czt - f_true;
        errors_improved_czt(j) = f_improved_czt - f_true;
        errors_rife(j) = f_rife - f_true;
        errors_mrife(j) = f_mrife - f_true;
        errors_irife(j) = f_irife - f_true;
        errors_iirife(j) = f_iirife - f_true;
    end

    % f. 计算当前SIR下的RMSE
    rmse_fft(i) = sqrt(mean(errors_fft.^2));
    rmse_czt(i) = sqrt(mean(errors_czt.^2));
    rmse_improved_czt(i) = sqrt(mean(errors_improved_czt.^2));
    rmse_rife(i) = sqrt(mean(errors_rife.^2));
    rmse_mrife(i) = sqrt(mean(errors_mrife.^2));
    rmse_irife(i) = sqrt(mean(errors_irife.^2));
    rmse_iirife(i) = sqrt(mean(errors_iirife.^2));
    
    % 计算标准差
    std_fft(i) = std(errors_fft);
    std_czt(i) = std(errors_czt);
    std_improved_czt(i) = std(errors_improved_czt);
    std_rife(i) = std(errors_rife);
    std_mrife(i) = std(errors_mrife);
    std_irife(i) = std(errors_irife);
    std_iirife(i) = std(errors_iirife);
end

fprintf('模拟完成。\n');

% --- 5. 绘制结果 ---
% 干扰类型显示字符串
interference_str = '';
switch interference_type
    case 'awgn'
        interference_str = '高斯白噪声';
    case 'sinusoidal'
        interference_str = sprintf('单频干扰(%.1f MHz)', f_interference/1e6);
    case 'chirp'
        interference_str = sprintf('线性调频干扰(%.1f-%.1f MHz)', f0_chirp/1e6, (f0_chirp + k_chirp*t(end))/1e6);
    case 'impulse'
        interference_str = sprintf('脉冲干扰(p=%.2f)', p_impulse);
end

% 在一个窗口中同时绘制RMSE和标准差图
figure('Position', [100, 100, 1200, 500]);

% 左侧：RMSE对比图
subplot(1, 2, 1);
semilogy(SIR_dB, rmse_fft, '-o', 'LineWidth', 1.5, 'DisplayName', 'FFT-Peak');
hold on;
semilogy(SIR_dB, rmse_czt, '-s', 'LineWidth', 1.5, 'DisplayName', 'CZT');
semilogy(SIR_dB, rmse_improved_czt, '-^', 'LineWidth', 1.5, 'DisplayName', '改进 CZT');
semilogy(SIR_dB, rmse_rife, '-d', 'LineWidth', 1.5, 'DisplayName', 'RIFE');
semilogy(SIR_dB, rmse_mrife, '-x', 'LineWidth', 1.5, 'DisplayName', 'MRIFE');
semilogy(SIR_dB, rmse_irife, '-+', 'LineWidth', 1.5, 'DisplayName', 'IRIFE');
semilogy(SIR_dB, rmse_iirife, '-*', 'LineWidth', 1.5, 'DisplayName', 'IIRIFE');
hold off;
grid on;
title('RMSE对比');
xlabel('SIR (dB)');
ylabel('RMSE (Hz)');
legend('show', 'Location', 'best');

% 右侧：标准差对比图
subplot(1, 2, 2);
semilogy(SIR_dB, std_fft, '-o', 'LineWidth', 1.5, 'DisplayName', 'FFT-Peak');
hold on;
semilogy(SIR_dB, std_czt, '-s', 'LineWidth', 1.5, 'DisplayName', 'CZT');
semilogy(SIR_dB, std_improved_czt, '-^', 'LineWidth', 1.5, 'DisplayName', '改进 CZT');
semilogy(SIR_dB, std_rife, '-d', 'LineWidth', 1.5, 'DisplayName', 'RIFE');
semilogy(SIR_dB, std_mrife, '-x', 'LineWidth', 1.5, 'DisplayName', 'MRIFE');
semilogy(SIR_dB, std_irife, '-+', 'LineWidth', 1.5, 'DisplayName', 'IRIFE');
semilogy(SIR_dB, std_iirife, '-*', 'LineWidth', 1.5, 'DisplayName', 'IIRIFE');
hold off;
grid on;
title('标准差对比');
xlabel('SIR (dB)');
ylabel('标准差 (Hz)');
legend('show', 'Location', 'best');

% 添加总标题
sgtitle(sprintf('七种频率估计算法性能对比 (频偏: %s, 干扰: %s, 窗函数: %s)', num2str(offset), interference_str, w_name));

