% =========================================================================
%   对比各算法在不同采样点数(N)下的运行时间。
% =========================================================================

clear;
close all;
clc;

% 添加 algorithms 目录到 MATLAB 路径
addpath('algorithms');

% --- 1. 仿真参数设置 ---
fs = 200e6; % 采样频率 (Hz)
N_values = 2.^(8:14); % 采样点数范围 (256, 512, ..., 16384)
num_n = length(N_values);

% 设置一个固定的频率和信噪比
f_center = 50e6; % 中心频率 (50 MHz)
offset = 0.4; % 相对频偏
SNR_dB = 20; % 固定高信噪比，以减少噪声对算法逻辑的影响

num_trials = 200; % 每个N值下的平均次数

% CZT 和 改进CZT 算法的参数
q = 1; % 频率细化区间大小控制参数
M = 64; % 频率细化倍数

% --- 2. 初始化结果存储变量 ---
time_fft = zeros(1, num_n);
time_czt = zeros(1, num_n);
time_improved_czt = zeros(1, num_n);
time_rife = zeros(1, num_n);
time_mrife = zeros(1, num_n);
time_irife = zeros(1, num_n);
time_iirife = zeros(1, num_n);

% --- 3. 执行计时模拟 ---
fprintf('开始运行时间对比模拟...\n');

for i = 1:num_n
    N = N_values(i);
    fprintf('正在处理 N = %d...\n', N);
    
    % a. 生成信号参数
    t = (0:N - 1) / fs;
    delta_f0 = fs / N;
    f_true = f_center + offset * delta_f0;
    
    % b. 生成纯净信号
    s_clean = exp(1j * 2 * pi * f_true * t);
    
    % c. 根据SNR计算噪声功率并生成噪声
    signal_power = mean(abs(s_clean) .^ 2);
    snr_linear = 10 ^ (SNR_dB / 10);
    noise_power = signal_power / snr_linear;
    noise = (randn(1, N) + 1j * randn(1, N)) * sqrt(noise_power / 2);
    
    % d. 生成带噪信号
    s_noisy = s_clean + noise;
    
    % e. 对每个算法重复计时取平均
    
    % FFT
    tic;
    for j = 1:num_trials, fft_est(s_noisy, fs); end
    time_fft(i) = toc / num_trials;
    
    % CZT
    tic;
    for j = 1:num_trials, czt_est(s_noisy, fs, q, M); end
    time_czt(i) = toc / num_trials;
    
    % 改进 CZT
    tic;
    for j = 1:num_trials, improved_czt_est(s_noisy, fs, q, M); end
    time_improved_czt(i) = toc / num_trials;
    
    % RIFE
    tic;
    for j = 1:num_trials, rife_est(s_noisy, fs); end
    time_rife(i) = toc / num_trials;
    
    % MRIFE
    tic;
    for j = 1:num_trials, mrife_est(s_noisy, fs); end
    time_mrife(i) = toc / num_trials;
    
    % IRIFE
    tic;
    for j = 1:num_trials, irife_est(s_noisy, fs); end
    time_irife(i) = toc / num_trials;
    
    % IIRIFE
    tic;
    for j = 1:num_trials, iirife_est(s_noisy, fs); end
    time_iirife(i) = toc / num_trials;
end

fprintf('模拟完成。\n');

% --- 4. 绘制结果 ---
figure;
semilogx(N_values, time_fft, '-o', 'LineWidth', 1.5, 'DisplayName', 'FFT-Peak');
hold on;
semilogx(N_values, time_czt, '-s', 'LineWidth', 1.5, 'DisplayName', 'CZT');
semilogx(N_values, time_improved_czt, '-^', 'LineWidth', 1.5, 'DisplayName', '改进 CZT');
semilogx(N_values, time_rife, '-d', 'LineWidth', 1.5, 'DisplayName', 'RIFE');
semilogx(N_values, time_mrife, '-x', 'LineWidth', 1.5, 'DisplayName', 'MRIFE');
semilogx(N_values, time_irife, '-+', 'LineWidth', 1.5, 'DisplayName', 'IRIFE');
semilogx(N_values, time_iirife, '-*', 'LineWidth', 1.5, 'DisplayName', 'IIRIFE');
hold off;
grid on;
title('七种频率估计算法运行时间对比');
xlabel('采样点数 N');
ylabel('平均运行时间 (秒)');
legend('show', 'Location', 'best');
