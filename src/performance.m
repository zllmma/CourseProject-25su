% =========================================================================
%               三种频率估计算法的性能对比脚本
% =========================================================================
%
% 目的:
%   对比直接FFT谱峰检测法、CZT法和改进CZT法在不同信噪比(SNR)下
%   的均方根误差(RMSE)和标准差(Std Dev)。
%
% 依赖:
%   - fft_est.m
%   - czt_est.m
%   - improved_czt_est.m
%
% =========================================================================

clear;
close all;
clc;

% --- 1. 仿真参数设置 ---
fs = 200e6;              % 采样频率 (Hz)
N = 1024;               % 采样点数
t = (0:N-1) / fs;       % 时间向量

% 设置一个非FFT整数倍的频率，以突显栅栏效应
f_true = 50.114e6;         % 真实信号频率 (Hz)

SNR_dB = -20:2:12;      % 信噪比范围 (dB)
num_trials = 1000;      % 每个SNR下的蒙特卡洛试验次数

% CZT 和 改进CZT 算法的参数
q = 1;                  % 频率细化区间大小控制参数
M = 64;                 % 频率细化倍数

% --- 2. 初始化结果存储变量 ---
num_snrs = length(SNR_dB);
rmse_fft = zeros(1, num_snrs);
rmse_czt = zeros(1, num_snrs);
rmse_improved_czt = zeros(1, num_snrs);
rmse_crlb = zeros(1, num_snrs);

std_fft = zeros(1, num_snrs);
std_czt = zeros(1, num_snrs);
std_improved_czt = zeros(1, num_snrs);

% --- 3. 执行蒙特卡洛模拟 ---
fprintf('开始蒙特卡洛模拟...\n');

% 并行执行每个SNR值的试验
% 使用 parfor 以加速处理
% 注意：parfor 需要 Parallel Computing Toolbox 支持
parfor i = 1:num_snrs
    snr_current_db = SNR_dB(i);
    fprintf('正在处理 SNR = %.0f dB...\n', snr_current_db);
    
    % 用于存储单次SNR下的所有试验误差
    errors_fft = zeros(1, num_trials);
    errors_czt = zeros(1, num_trials);
    errors_improved_czt = zeros(1, num_trials);
    phases = 2 * pi * randn(1, num_trials); % 随机相位用于每次试验

    for j = 1:num_trials
        % a. 生成纯净信号 (使用复正弦信号)
        % 为增加随机性，每次试验使用不同的初始相位
        phi = phases(j);
        s_clean = exp(1j * (2 * pi * f_true * t + phi));

        % b. 根据SNR计算噪声功率并生成噪声
        signal_power = mean(abs(s_clean).^2);
        snr_linear = 10^(snr_current_db / 10);
        noise_power = signal_power / snr_linear;
        
        % 生成复高斯白噪声
        noise = (randn(1, N) + 1j * randn(1, N)) * sqrt(noise_power / 2);

        % c. 生成带噪信号
        s_noisy = s_clean + noise;

        % d. 使用三种算法进行频率估计
        f_fft = fft_est(s_noisy, fs);
        f_czt = czt_est(s_noisy, fs, q, M);
        f_improved_czt = improved_czt_est(s_noisy, fs, q, M);

        % e. 计算并存储误差
        errors_fft(j) = f_fft - f_true;
        errors_czt(j) = f_czt - f_true;
        errors_improved_czt(j) = f_improved_czt - f_true;
    end

    % f. 计算当前SNR下的RMSE和标准差
    rmse_fft(i) = sqrt(mean(errors_fft.^2));
    rmse_czt(i) = sqrt(mean(errors_czt.^2));
    rmse_improved_czt(i) = sqrt(mean(errors_improved_czt.^2));
    snr_linear = 10^(snr_current_db / 10);
    rmse_crlb(i) = sqrt((6 * fs^2) / ((2*pi)^2 * snr_linear * N * (N^2 - 1)));

    std_fft(i) = std(errors_fft);
    std_czt(i) = std(errors_czt);
    std_improved_czt(i) = std(errors_improved_czt);
end

fprintf('模拟完成。\n');

% --- 4. 绘制结果 ---
% 绘制RMSE对比图
figure;
semilogy(SNR_dB, rmse_fft, '-o', 'LineWidth', 1.5, 'DisplayName', 'FFT-Peak');
hold on;
semilogy(SNR_dB, rmse_czt, '-s', 'LineWidth', 1.5, 'DisplayName', 'CZT');
semilogy(SNR_dB, rmse_improved_czt, '-^', 'LineWidth', 1.5, 'DisplayName', '改进 CZT');
semilogy(SNR_dB, rmse_crlb, 'k--', 'LineWidth', 2, 'DisplayName', 'CRLB');
hold off;
grid on;
title('三种频率估计算法的RMSE对比');
xlabel('SNR (dB)');
ylabel('RMSE (Hz)');
legend('show', 'Location', 'northeast');

% 绘制标准差对比图
figure;
semilogy(SNR_dB, std_fft, '-o', 'LineWidth', 1.5, 'DisplayName', 'FFT');
hold on;
semilogy(SNR_dB, std_czt, '-s', 'LineWidth', 1.5, 'DisplayName', 'CZT');
semilogy(SNR_dB, std_improved_czt, '-^', 'LineWidth', 1.5, 'DisplayName', '改进 CZT');
hold off;
grid on;
title('三种频率估计算法的标准差对比');
xlabel('SNR (dB)');
ylabel('SD (Hz)');
legend('show', 'Location', 'northeast');