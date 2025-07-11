% =========================================================================
%      在固定低信噪比下，RMSE随相对频偏变化的对比脚本 (支持窗函数)
% =========================================================================
%
% 新增功能:
%   1. 可选择窗函数：'none'(不使用窗), 'rect'(矩形窗), 'hann'(Hann窗), 
%      'hamming'(Hamming窗), 'blackman'(Blackman窗), 'kaiser'(Kaiser窗)
%   2. CRLB仅在无窗函数时显示（因为CRLB理论适用于无窗情况）
%
% 使用示例:
%   window_type = 'hann';  % 窗函数类型
%   SNR_dB = 10;           % 信噪比 (dB)
% =========================================================================

clear;
close all;
clc;

% 添加 algorithms 目录到 MATLAB 路径
addpath('algorithms');

% --- 0. 用户参数选择 (修改此处) ---
window_type = 'hamming';              % 可选：'none', 'rect', 'hann', 'hamming', 'blackman', 'kaiser'

% --- 1. 模拟参数设置 ---
fs = 200e6; % 采样频率 (200 MHz)
N = 1024; % 采样点数
t = (0:N - 1) / fs; % 时间向量
A = 1.0; % 信号幅度

f_center = 50e6; % 中心频率 (50 MHz)
delta_f0 = fs / N; % 频率分辨率

SNR_dB = 0; % 固定信噪比 (dB) - 更改为更合理的值

% 设置相对频偏的范围
relative_offsets = -0.5:0.05:0.5;

num_trials = 1000; % 每个频偏下的蒙特卡洛试验次数

% CZT 参数
q = 1;
M = 64;

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
        error('未知窗类型: %s。支持的窗类型: none, rect, hann, hamming, blackman, kaiser', window_type);
end

% 窗函数归一化 (保持信号功率不变)
% 对于'none'选项不进行归一化，保持原始信号的幅度
if ~strcmp(window_type, 'none')
    w = w / sqrt(mean(w.^2));  % 功率归一化：保持信号功率不变
end

% --- 3. 初始化结果存储变量 ---
num_offsets = length(relative_offsets);
rmse_fft = zeros(1, num_offsets);
rmse_czt = zeros(1, num_offsets);
rmse_improved_czt = zeros(1, num_offsets);
rmse_rife = zeros(1, num_offsets);
rmse_mrife = zeros(1, num_offsets);
rmse_irife = zeros(1, num_offsets);
rmse_iirife = zeros(1, num_offsets);
rmse_crlb = zeros(1, num_offsets); % 添加CRLB存储变量

% 添加标准差存储变量
std_fft = zeros(1, num_offsets);
std_czt = zeros(1, num_offsets);
std_improved_czt = zeros(1, num_offsets);
std_rife = zeros(1, num_offsets);
std_mrife = zeros(1, num_offsets);
std_irife = zeros(1, num_offsets);
std_iirife = zeros(1, num_offsets);

% --- 4. 执行主循环 (遍历所有频偏) ---
fprintf('开始在 SNR = %.0f dB 下进行模拟，窗函数: %s...\n', SNR_dB, w_name);

% 提前计算噪声参数
snr_linear = 10 ^ (SNR_dB / 10);
signal_power = A ^ 2; % 复信号功率 (假设无窗或窗函数已归一化)
noise_power = signal_power / snr_linear;
noise_std_per_component = sqrt(noise_power / 2);

% 计算CRLB (Cramér-Rao Lower Bound) - 仅适用于无窗情况
% CRLB = sqrt((6 * fs^2) / ((2*pi)^2 * snr_linear * N * (N^2 - 1)))
% CRLB与频偏无关，仅与SNR、采样点数和采样频率有关
if strcmp(window_type, 'none')
    rmse_crlb(:) = sqrt((6 * fs ^ 2) / ((2 * pi) ^ 2 * snr_linear * N * (N ^ 2 - 1)));
    show_crlb = true;
else
    show_crlb = false;
end

% 使用 parfor 并行处理每个相对频偏
% 以加速模拟过程
% 注意：parfor 需要 Parallel Computing Toolbox 支持
parfor i = 1:num_offsets
    current_offset = relative_offsets(i);

    % 根据当前频偏计算真实频率
    f_true = f_center + current_offset * delta_f0;

    fprintf('正在处理相对频偏: %.2f, 真实频率: %.4f MHz\n', current_offset, f_true / 1e6);

    % 用于存储单次频偏下的所有试验误差
    errors_fft_mc = zeros(1, num_trials);
    errors_czt_mc = zeros(1, num_trials);
    errors_improved_czt_mc = zeros(1, num_trials);
    errors_rife_mc = zeros(1, num_trials);
    errors_mrife_mc = zeros(1, num_trials);
    errors_irife_mc = zeros(1, num_trials);
    errors_iirife_mc = zeros(1, num_trials);
    phases = 2 * pi * rand(1, num_trials); % 随机相位用于每次试验

    % 蒙特卡洛模拟
    for j = 1:num_trials
        % a. 生成纯净信号
        phi = phases(j); % 随机相位
        s_clean = A * exp(1j * (2 * pi * f_true * t + phi));
        
        % b. 先对纯净信号应用窗函数 (用于正确计算信号功率)
        s_windowed_clean = s_clean .* w;
        
        % c. 基于加窗后的信号功率生成噪声
        windowed_signal_power = mean(abs(s_windowed_clean).^2);
        actual_snr_linear = 10 ^ (SNR_dB / 10);
        actual_noise_power = windowed_signal_power / actual_snr_linear;
        actual_noise_std = sqrt(actual_noise_power / 2);
        
        % d. 生成噪声并应用相同的窗函数
        noise = (randn(1, N) + 1j * randn(1, N)) * actual_noise_std;
        noise_windowed = noise .* w;
        
        % e. 最终的含噪信号
        s_noisy = s_windowed_clean + noise_windowed;

        % f. 使用七种算法进行频率估计
        f_fft = fft_est(s_noisy, fs);
        f_czt = czt_est(s_noisy, fs, q, M);
        f_improved_czt = improved_czt_est(s_noisy, fs, q, M);
        f_rife = rife_est(s_noisy, fs);
        f_mrife = mrife_est(s_noisy, fs);
        f_irife = irife_est(s_noisy, fs);
        f_iirife = iirife_est(s_noisy, fs);

        % g. 计算并存储误差
        errors_fft_mc(j) = f_fft - f_true;
        errors_czt_mc(j) = f_czt - f_true;
        errors_improved_czt_mc(j) = f_improved_czt - f_true;
        errors_rife_mc(j) = f_rife - f_true;
        errors_mrife_mc(j) = f_mrife - f_true;
        errors_irife_mc(j) = f_irife - f_true;
        errors_iirife_mc(j) = f_iirife - f_true;
    end

    % h. 计算当前频偏下的RMSE
    rmse_fft(i) = sqrt(mean(errors_fft_mc .^ 2));
    rmse_czt(i) = sqrt(mean(errors_czt_mc .^ 2));
    rmse_improved_czt(i) = sqrt(mean(errors_improved_czt_mc .^ 2));
    rmse_rife(i) = sqrt(mean(errors_rife_mc .^ 2));
    rmse_mrife(i) = sqrt(mean(errors_mrife_mc .^ 2));
    rmse_irife(i) = sqrt(mean(errors_irife_mc .^ 2));
    rmse_iirife(i) = sqrt(mean(errors_iirife_mc .^ 2));
    
    % 计算标准差
    std_fft(i) = std(errors_fft_mc);
    std_czt(i) = std(errors_czt_mc);
    std_improved_czt(i) = std(errors_improved_czt_mc);
    std_rife(i) = std(errors_rife_mc);
    std_mrife(i) = std(errors_mrife_mc);
    std_irife(i) = std(errors_irife_mc);
    std_iirife(i) = std(errors_iirife_mc);
end

fprintf('模拟完成。\n');

% --- 4. 在一个窗口中绘制RMSE和标准差图 ---
figure('Position', [100, 100, 1200, 500]);

% 左侧：RMSE对比图
subplot(1, 2, 1);
semilogy(relative_offsets, rmse_fft, '-o', 'LineWidth', 1.5, 'DisplayName', 'FFT-Peak');
hold on;
semilogy(relative_offsets, rmse_czt, '-s', 'LineWidth', 1.5, 'DisplayName', 'CZT');
semilogy(relative_offsets, rmse_improved_czt, '-^', 'LineWidth', 1.5, 'DisplayName', '改进 CZT');
semilogy(relative_offsets, rmse_rife, '-d', 'LineWidth', 1.5, 'DisplayName', 'Rife');
semilogy(relative_offsets, rmse_mrife, '-x', 'LineWidth', 1.5, 'DisplayName', 'M-Rife');
semilogy(relative_offsets, rmse_irife, '-+', 'LineWidth', 1.5, 'DisplayName', 'I-Rife');
semilogy(relative_offsets, rmse_iirife, '-*', 'LineWidth', 1.5, 'DisplayName', 'IIRife');

% 只在无窗函数时显示CRLB理论下界
if show_crlb
    semilogy(relative_offsets, rmse_crlb, 'k--', 'LineWidth', 2, 'DisplayName', 'CRLB', 'HandleVisibility', 'on');
end

hold off;
grid on;
title(['SNR = ' num2str(SNR_dB) ' dB 时, RMSE随相对频偏的变化']);
xlabel('相对频偏 (单位: 频率分辨率)');
ylabel('均方根误差 (RMSE, Hz)');
legend('show', 'Location', 'best');

% 右侧：标准差对比图
subplot(1, 2, 2);
semilogy(relative_offsets, std_fft, '-o', 'LineWidth', 1.5, 'DisplayName', 'FFT-Peak');
hold on;
semilogy(relative_offsets, std_czt, '-s', 'LineWidth', 1.5, 'DisplayName', 'CZT');
semilogy(relative_offsets, std_improved_czt, '-^', 'LineWidth', 1.5, 'DisplayName', '改进 CZT');
semilogy(relative_offsets, std_rife, '-d', 'LineWidth', 1.5, 'DisplayName', 'Rife');
semilogy(relative_offsets, std_mrife, '-x', 'LineWidth', 1.5, 'DisplayName', 'M-Rife');
semilogy(relative_offsets, std_irife, '-+', 'LineWidth', 1.5, 'DisplayName', 'I-Rife');
semilogy(relative_offsets, std_iirife, '-*', 'LineWidth', 1.5, 'DisplayName', 'IIRife');
hold off;
grid on;
title(['SNR = ' num2str(SNR_dB) ' dB 时, 标准差随相对频偏的变化']);
xlabel('相对频偏 (单位: 频率分辨率)');
ylabel('标准差 (Hz)');
legend('show', 'Location', 'best');

% 调整整个窗口布局
sgtitle(sprintf('七种频率估计算法性能对比 (窗函数: %s, 噪声: AWGN)', w_name));
