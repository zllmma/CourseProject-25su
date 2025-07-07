% =========================================================================
%      在固定低信噪比下，RMSE随相对频偏变化的对比脚本
% =========================================================================
clear;
close all;
clc;

% --- 1. 模拟参数设置 ---
fs = 200e6;             % 采样频率 (200 MHz)
N = 1024;               % 采样点数
t = (0:N-1) / fs;       % 时间向量
A = 1.0;                % 信号幅度

f_center = 50e6;        % 中心频率 (50 MHz)
delta_f0 = fs / N;      % 频率分辨率

SNR_dB = -8;            % 固定信噪比 (dB)

% 设置相对频偏的范围
relative_offsets = -0.5:0.1:0.5;

num_trials = 1000;      % 每个频偏下的蒙特卡洛试验次数

% CZT 参数
q = 1;
M = 64;

% --- 2. 初始化结果存储变量 ---
num_offsets = length(relative_offsets);
rmse_fft = zeros(1, num_offsets);
rmse_czt = zeros(1, num_offsets);
rmse_improved_czt = zeros(1, num_offsets);

% --- 3. 执行主循环 (遍历所有频偏) ---
fprintf('开始在 SNR = %.0f dB 下进行模拟...\n', SNR_dB);

% 提前计算噪声参数
snr_linear = 10^(SNR_dB / 10);
signal_power = A^2; % 复信号功率
noise_power = signal_power / snr_linear;
noise_std_per_component = sqrt(noise_power / 2);

parfor i = 1:num_offsets
    current_offset = relative_offsets(i);
    
    % 根据当前频偏计算真实频率
    f_true = f_center + current_offset * delta_f0;
    
    fprintf('正在处理相对频偏: %.1f, 真实频率: %.4f MHz\n', current_offset, f_true/1e6);
    
    % 用于存储单次频偏下的所有试验误差
    errors_fft_mc = zeros(1, num_trials);
    errors_czt_mc = zeros(1, num_trials);
    errors_improved_czt_mc = zeros(1, num_trials);
    phases = 2 * pi * rand(1, num_trials); % 随机相位用于每次试验
    
    % Monte Carlo 模拟
    for j = 1:num_trials
        % a. 生成信号和噪声
        phi = phases(j); % 随机相位
        s_clean = A * exp(1j * (2 * pi * f_true * t + phi));
        noise = (randn(1, N) + 1j * randn(1, N)) * noise_std_per_component;
        s_noisy = s_clean + noise;

        % b. 使用三种算法进行频率估计
        f_fft = fft_est(s_noisy, fs);
        f_czt = czt_est(s_noisy, fs, q, M);
        f_improved_czt = improved_czt_est(s_noisy, fs, q, M);

        % c. 计算并存储误差
        errors_fft_mc(j) = f_fft - f_true;
        errors_czt_mc(j) = f_czt - f_true;
        errors_improved_czt_mc(j) = f_improved_czt - f_true;
    end
    
    % d. 计算当前频偏下的RMSE
    rmse_fft(i) = sqrt(mean(errors_fft_mc.^2));
    rmse_czt(i) = sqrt(mean(errors_czt_mc.^2));
    rmse_improved_czt(i) = sqrt(mean(errors_improved_czt_mc.^2));
end

fprintf('模拟完成。\n');

% --- 4. 绘制结果 ---
figure;
semilogy(relative_offsets, rmse_fft, '-o', 'LineWidth', 1.5, 'DisplayName', 'FFT-Peak');
hold on;
semilogy(relative_offsets, rmse_czt, '-s', 'LineWidth', 1.5, 'DisplayName', 'CZT');
semilogy(relative_offsets, rmse_improved_czt, '-^', 'LineWidth', 1.5, 'DisplayName', '改进 CZT');
hold off;
grid on;
title(['SNR = ' num2str(SNR_dB) ' dB 时, RMSE随相对频偏的变化']);
xlabel('相对频偏 (单位: 频率分辨率)');
ylabel('均方根误差 (RMSE, Hz)');
legend('show', 'Location', 'best');
set(gca, 'FontSize', 12);