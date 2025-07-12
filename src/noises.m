% =========================================================================
%               七种频率估计算法的性能对比脚本 (支持多种噪声和窗函数)
% =========================================================================
%
% 功能:
%   1. 可选择噪声类型：
%      - 'awgn'(高斯白噪声)
%      - 'uniform'(均匀分布噪声)
%      - 'impulse'(脉冲噪声)
%      - 'rayleigh'(瑞利噪声)
%      - 'poisson'(泊松噪声)
%   2. 可选择窗函数：'none'(不使用窗), 'hann'(Hann窗), 'hamming'(Hamming窗), 'blackman'(Blackman窗), 'kaiser'(Kaiser窗)
%
% 使用示例:
%   noise_type = 'awgn';      % 高斯白噪声 (默认)
%   window_type = 'hann';     % 可选：'none', 'hann', 'hamming', 'blackman', 'kaiser'
%   p_impulse = 0.1;          % 脉冲噪声出现概率
%   poisson_lambda = 10;      % 泊松噪声参数
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

% --- 0. 用户参数选择 (修改此处) ---
noise_type = 'awgn';     % 可选：'awgn', 'uniform', 'impulse', 'rayleigh', 'poisson'
window_type = 'hann';    % 可选：'none', 'rect', 'hann', 'hamming', 'blackman'
p_impulse = 0.1;         % 脉冲噪声出现概率 (仅noise_type='impulse'时有效)
poisson_lambda = 10;     % 泊松分布参数 (仅noise_type='poisson'时有效)

% --- 1. 仿真参数设置 ---
fs = 200e6;              % 采样频率 (Hz)
N = 1024;                % 采样点数
t = (0:N-1) / fs;        % 时间向量

% 设置一个非FFT整数倍的频率，以突显栅栏效应
f_center = 50e6; % 中心频率 (50 MHz)
delta_f0 = fs / N; % 频率分辨率 (Hz)
offset = 0.1; % 相对频偏
f_true = f_center + offset * delta_f0; % 真实信号频率 (Hz)

SNR_dB = -20:2:20;       % 信噪比范围 (dB)
num_trials = 1000;       % 每个SNR下的蒙特卡洛试验次数

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

% --- 3. 初始化结果存储变量 ---
num_snrs = length(SNR_dB);
rmse_fft = zeros(1, num_snrs);
rmse_czt = zeros(1, num_snrs);
rmse_improved_czt = zeros(1, num_snrs);
rmse_rife = zeros(1, num_snrs);
rmse_mrife = zeros(1, num_snrs);
rmse_irife = zeros(1, num_snrs);
rmse_iirife = zeros(1, num_snrs);

% 添加标准差存储变量
std_fft = zeros(1, num_snrs);
std_czt = zeros(1, num_snrs);
std_improved_czt = zeros(1, num_snrs);
std_rife = zeros(1, num_snrs);
std_mrife = zeros(1, num_snrs);
std_irife = zeros(1, num_snrs);
std_iirife = zeros(1, num_snrs);

% --- 4. 执行蒙特卡洛模拟 ---
fprintf('开始蒙特卡洛模拟，噪声类型: %s, 窗类型: %s...\n', noise_type, w_name);

% 并行执行每个SNR值的试验
parfor i = 1:num_snrs
    snr_current_db = SNR_dB(i);
    fprintf('正在处理 SNR = %.0f dB...\n', snr_current_db);
    
    % 用于存储单次SNR下的所有试验误差
    errors_fft = zeros(1, num_trials);
    errors_czt = zeros(1, num_trials);
    errors_improved_czt = zeros(1, num_trials);
    errors_rife = zeros(1, num_trials);
    errors_mrife = zeros(1, num_trials);
    errors_irife = zeros(1, num_trials);
    errors_iirife = zeros(1, num_trials);
    phases = 2 * pi * randn(1, num_trials); % 随机相位用于每次试验
    
    for j = 1:num_trials
        % a. 生成纯净信号 (使用复正弦信号)
        phi = phases(j);
        s_clean = exp(1j * (2 * pi * f_true * t + phi));
        
        % b. 先对纯净信号加窗，然后基于加窗后的信号功率计算噪声
        s_windowed_clean = s_clean .* w;
        windowed_signal_power = mean(abs(s_windowed_clean).^2);
        snr_linear = 10^(snr_current_db / 10);
        noise_power = windowed_signal_power / snr_linear;
        
        % 根据噪声类型生成噪声
        % 初始化噪声向量
        noise_real = zeros(1, N);
        noise_imag = zeros(1, N);
        
        switch noise_type
            case 'awgn'
                noise_real = randn(1, N) * sqrt(noise_power/2);
                noise_imag = randn(1, N) * sqrt(noise_power/2);
            case 'uniform'
                a = sqrt(3*noise_power);
                noise_real = (rand(1, N)*2 - 1) * a;
                noise_imag = (rand(1, N)*2 - 1) * a;
            case 'impulse'
                sigma_impulse = sqrt(noise_power/(2*p_impulse));
                impulse_real = double(rand(1, N) < p_impulse) .* ...
                    (sigma_impulse * randn(1, N));
                impulse_imag = double(rand(1, N) < p_impulse) .* ...
                    (sigma_impulse * randn(1, N));
                noise_real = impulse_real;
                noise_imag = impulse_imag;
            case 'rayleigh'
                % 瑞利噪声
                % 调整参数使得噪声功率为noise_power
                sigma = sqrt(noise_power/2);
                % 生成瑞利分布的幅度
                mag = sigma * sqrt(-2*log(rand(1, N)));
                % 均匀分布的相位
                phase = 2*pi*rand(1, N);
                % 转换为实部和虚部
                noise_real = mag .* cos(phase);
                noise_imag = mag .* sin(phase);
            case 'poisson'
                % 泊松噪声 - 调整信号强度使总功率符合SNR要求
                % 泊松过程往往用于建模光子计数过程
                scaling = sqrt(noise_power/poisson_lambda);
                % 生成泊松随机变量
                noise_real = (poissrnd(poisson_lambda, 1, N) - poisson_lambda) * scaling;
                noise_imag = (poissrnd(poisson_lambda, 1, N) - poisson_lambda) * scaling;
        end

        noise = noise_real + 1j * noise_imag;

        % c. 合成最终的带噪信号（纯净信号已经加窗）
        s_windowed = s_windowed_clean + noise;
        
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
    
    % f. 计算当前SNR下的RMSE
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
% 噪声类型显示字符串
noise_str = '';
switch noise_type
    case 'awgn'
        noise_str = '高斯白噪声';
    case 'uniform'
        noise_str = '均匀分布噪声';
    case 'impulse'
        noise_str = sprintf('脉冲噪声(p=%.2f)', p_impulse);
    case 'rayleigh'
        noise_str = '瑞利噪声';
    case 'poisson'
        noise_str = sprintf('泊松噪声(λ=%.1f)', poisson_lambda);
end

% 在一个窗口中同时绘制RMSE和标准差图
figure('Position', [100, 100, 1200, 500]);

% 左侧：RMSE对比图
subplot(1, 2, 1);
semilogy(SNR_dB, rmse_fft, '-o', 'LineWidth', 1.5, 'DisplayName', 'FFT-Peak');
hold on;
semilogy(SNR_dB, rmse_czt, '-s', 'LineWidth', 1.5, 'DisplayName', 'CZT');
semilogy(SNR_dB, rmse_improved_czt, '-^', 'LineWidth', 1.5, 'DisplayName', '改进 CZT');
semilogy(SNR_dB, rmse_rife, '-d', 'LineWidth', 1.5, 'DisplayName', 'RIFE');
semilogy(SNR_dB, rmse_mrife, '-x', 'LineWidth', 1.5, 'DisplayName', 'MRIFE');
semilogy(SNR_dB, rmse_irife, '-+', 'LineWidth', 1.5, 'DisplayName', 'IRIFE');
semilogy(SNR_dB, rmse_iirife, '-*', 'LineWidth', 1.5, 'DisplayName', 'IIRIFE');

% 添加CRLB理论下界曲线 (仅当噪声类型为AWGN且不使用窗函数时)
if strcmp(noise_type, 'awgn') && strcmp(window_type, 'none')
    % 计算CRLB曲线
    snr_linear = 10.^(SNR_dB / 10);
    % 单频信号的CRLB
    crlb = sqrt(3 / (2 * pi^2 * N * (N^2 - 1))) * fs ./ sqrt(snr_linear);
    % 在图上绘制CRLB曲线
    semilogy(SNR_dB, crlb, '--k', 'LineWidth', 2, 'DisplayName', 'CRLB');
end

hold off;
grid on;
title('RMSE对比');
xlabel('SNR (dB)');
ylabel('RMSE (Hz)');
legend('show', 'Location', 'best');

% 右侧：标准差对比图
subplot(1, 2, 2);
semilogy(SNR_dB, std_fft, '-o', 'LineWidth', 1.5, 'DisplayName', 'FFT-Peak');
hold on;
semilogy(SNR_dB, std_czt, '-s', 'LineWidth', 1.5, 'DisplayName', 'CZT');
semilogy(SNR_dB, std_improved_czt, '-^', 'LineWidth', 1.5, 'DisplayName', '改进 CZT');
semilogy(SNR_dB, std_rife, '-d', 'LineWidth', 1.5, 'DisplayName', 'RIFE');
semilogy(SNR_dB, std_mrife, '-x', 'LineWidth', 1.5, 'DisplayName', 'MRIFE');
semilogy(SNR_dB, std_irife, '-+', 'LineWidth', 1.5, 'DisplayName', 'IRIFE');
semilogy(SNR_dB, std_iirife, '-*', 'LineWidth', 1.5, 'DisplayName', 'IIRIFE');
hold off;
grid on;
title('标准差对比');
xlabel('SNR (dB)');
ylabel('标准差 (Hz)');
legend('show', 'Location', 'best');

% 添加总标题
sgtitle(sprintf('七种频率估计算法性能对比 (频偏: %s, 噪声: %s, 窗函数: %s)', num2str(offset), noise_str, w_name));