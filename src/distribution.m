% =========================================================================
%      在固定信噪比下，七种频率估计算法的分布对比脚本
% =========================================================================
clear;
close all;
clc;

% 添加 algorithms 目录到 MATLAB 路径
% 首先检查路径是否存在，然后添加
if exist('src/algorithms', 'dir')
    addpath('src/algorithms');
elseif exist('algorithms', 'dir')
    addpath('algorithms');
else
    error('无法找到算法目录，请确保 algorithms 目录存在');
end

% --- 1. 模拟参数设置 ---
fs = 200e6; % 采样频率 (Hz)
N = 1024; % 采样点数
t = (0:N - 1) / fs; % 时间向量
A = 1;
f_center = 50e6; % 中心频率 (50 MHz)
delta_f0 = fs / N; % 频率分辨率 (Hz)
offset = 0.1; % 相对频偏
f_true = f_center + offset * delta_f0; % 真实信号频率 (Hz)

% 【关键】选择一个固定的信噪比进行分析
SNR_fixed_dB = -8; % 固定信噪比 (dB)

num_trials = 10000; % 大量重复试验以获得平滑的分布

% CZT 和 改进CZT 算法的参数
q = 1; % 频率细化区间大小控制参数
M = 64; % 频率细化倍数

% --- 2. 初始化结果存储变量 ---
estimates_fft = zeros(1, num_trials);
estimates_czt = zeros(1, num_trials);
estimates_improved_czt = zeros(1, num_trials);
estimates_rife = zeros(1, num_trials);
estimates_mrife = zeros(1, num_trials);
estimates_irife = zeros(1, num_trials);
estimates_iirife = zeros(1, num_trials);

% 计算一次噪声功率即可
snr_linear = 10 ^ (SNR_fixed_dB / 10);
signal_power = A ^ 2;
noise_power = signal_power / snr_linear;
noise_std = sqrt(noise_power / 2);
phases = 2 * pi * rand(1, num_trials); % 随机相位用于每次试验

% 并行执行蒙特卡洛模拟
% 使用 parfor 以加速处理
% 注意：parfor 需要 Parallel Computing Toolbox 支持
parfor i = 1:num_trials
    % a. 生成纯净信号 (每次相位随机)
    phi = phases(i);
    s_clean = exp(1j * (2 * pi * f_true * t + phi));

    % b. 生成噪声
    noise = (randn(1, N) + 1j * randn(1, N)) * noise_std;

    % c. 生成带噪信号
    s_noisy = s_clean + noise;

    % d. 使用七种算法进行频率估计
    estimates_fft(i) = fft_est(s_noisy, fs);
    estimates_czt(i) = czt_est(s_noisy, fs, q, M);
    estimates_improved_czt(i) = improved_czt_est(s_noisy, fs, q, M);
    estimates_rife(i) = rife_est(s_noisy, fs);
    estimates_mrife(i) = mrife_est(s_noisy, fs);
    estimates_irife(i) = irife_est(s_noisy, fs);
    estimates_iirife(i) = iirife_est(s_noisy, fs);

    % 打印进度
    if mod(i, 1000) == 0
        fprintf('已完成 %d / %d 次试验...\n', i, num_trials);
    end

end

fprintf('模拟完成。\n');

% --- 4. 绘制分布直方图 (修改为三个独立的子图)

% 为了方便比较，统一所有子图的X轴范围
min_freq = min([estimates_fft, estimates_czt, estimates_improved_czt, estimates_rife, estimates_mrife, estimates_irife, estimates_iirife]);
max_freq = max([estimates_fft, estimates_czt, estimates_improved_czt, estimates_rife, estimates_mrife, estimates_irife, estimates_iirife]);
x_limits = [min_freq - 1e3, max_freq + 1e3]; % 稍微留出一些边距

% 第一个子图: FFT 方法
subplot(7, 1, 1);
histogram(estimates_fft, 50, 'Normalization', 'pdf', 'FaceColor', 'b', 'FaceAlpha', 0.75);
hold on;
xline(f_true, 'r--', 'LineWidth', 2);
hold off;
title('FFT 方法');
ylabel('概率密度');
grid on;
xlim(x_limits); % 应用统一的X轴范围

% 第二个子图: CZT 方法
subplot(7, 1, 2);
histogram(estimates_czt, 50, 'Normalization', 'pdf', 'FaceColor', 'g', 'FaceAlpha', 0.75);
hold on;
xline(f_true, 'r--', 'LineWidth', 2);
hold off;
title('CZT 方法');
ylabel('概率密度');
grid on;
xlim(x_limits); % 应用统一的X轴范围

% 第三个子图: 改进 CZT 方法
subplot(7, 1, 3);
histogram(estimates_improved_czt, 50, 'Normalization', 'pdf', 'FaceColor', 'm', 'FaceAlpha', 0.75);
hold on;
xline(f_true, 'r--', 'LineWidth', 2, 'DisplayName', '真实频率');
hold off;
title('改进 CZT 方法');
xlabel('估计的频率 (Hz)');
ylabel('概率密度');
grid on;
xlim(x_limits); % 应用统一的X轴范围

% 第四个子图: RIFE 方法
subplot(7, 1, 4);
histogram(estimates_rife, 50, 'Normalization', 'pdf', 'FaceColor', 'c', 'FaceAlpha', 0.75);
hold on;
xline(f_true, 'r--', 'LineWidth', 2);
hold off;
title('RIFE 方法');
ylabel('概率密度');
grid on;
xlim(x_limits); % 应用统一的X轴范围

% 第五个子图: MRIFE 方法
subplot(7, 1, 5);
histogram(estimates_mrife, 50, 'Normalization', 'pdf', 'FaceColor', 'y', 'FaceAlpha', 0.75);
hold on;
xline(f_true, 'r--', 'LineWidth', 2);
hold off;
title('MRIFE 方法');
ylabel('概率密度');
grid on;
xlim(x_limits); % 应用统一的X轴范围

% 第六个子图: IRIFE 方法
subplot(7, 1, 6);
histogram(estimates_irife, 50, 'Normalization', 'pdf', 'FaceColor', 'k', 'FaceAlpha', 0.75);
hold on;
xline(f_true, 'r--', 'LineWidth', 2);
hold off;
title('IRIFE 方法');
ylabel('概率密度');
grid on;
xlim(x_limits); % 应用统一的X轴范围

% 第七个子图: IIRIFE 方法
subplot(7, 1, 7);
histogram(estimates_iirife, 50, 'Normalization', 'pdf', 'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.75);
hold on;
xline(f_true, 'r--', 'LineWidth', 2);
hold off;
title('IIRIFE 方法');
ylabel('概率密度');
grid on;
xlim(x_limits); % 应用统一的X轴范围

% 添加一个总标题
sgtitle(['SNR = ' num2str(SNR_fixed_dB) ' dB 时, 各算法频率估计分布'], 'FontSize', 14, 'FontWeight', 'bold');
