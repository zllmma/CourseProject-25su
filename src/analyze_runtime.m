% 运行时间分析脚本
clear; clc;

% 添加 algorithms 目录到 MATLAB 路径
addpath('algorithms');

% 参数设置
fs = 200e6;
N_values = 2.^(8:14); % 256, 512, 1024, 2048, 4096, 8192, 16384
num_n = length(N_values);

f_center = 50e6;
offset = 0.4;
SNR_dB = 20;
num_trials = 50; % 减少试验次数以快速获得结果

% CZT 参数
q = 1;
M = 64;

% 初始化结果存储
time_fft = zeros(1, num_n);
time_czt = zeros(1, num_n);
time_improved_czt = zeros(1, num_n);
time_rife = zeros(1, num_n);
time_mrife = zeros(1, num_n);
time_irife = zeros(1, num_n);
time_iirife = zeros(1, num_n);

fprintf('=== 频率估计算法运行时间对比分析 ===\n');
fprintf('采样频率: %.0f MHz, SNR: %d dB, 平均次数: %d\n\n', fs/1e6, SNR_dB, num_trials);

for i = 1:num_n
    N = N_values(i);
    fprintf('正在测试 N = %d...', N);
    
    % 生成测试信号
    t = (0:N-1) / fs;
    delta_f0 = fs / N;
    f_true = f_center + offset * delta_f0;
    
    s_clean = exp(1j * 2 * pi * f_true * t);
    signal_power = mean(abs(s_clean).^2);
    snr_linear = 10^(SNR_dB/10);
    noise_power = signal_power / snr_linear;
    noise = (randn(1, N) + 1j * randn(1, N)) * sqrt(noise_power/2);
    s_noisy = s_clean + noise;
    
    % 计时测试
    % FFT
    tic;
    for j = 1:num_trials, fft_est(s_noisy, fs); end
    time_fft(i) = toc / num_trials;
    
    % CZT
    tic;
    for j = 1:num_trials, czt_est(s_noisy, fs, q, M); end
    time_czt(i) = toc / num_trials;
    
    % 改进CZT
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
    
    fprintf(' 完成\n');
end

% 显示结果表格
fprintf('\n=== 运行时间结果 (毫秒) ===\n');
fprintf('%8s %8s %8s %8s %8s %8s %8s %8s\n', 'N', 'FFT', 'CZT', '改进CZT', 'RIFE', 'MRIFE', 'IRIFE', 'IIRIFE');
fprintf('%s\n', repmat('-', 1, 80));

for i = 1:num_n
    fprintf('%8d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
        N_values(i), ...
        time_fft(i)*1000, time_czt(i)*1000, time_improved_czt(i)*1000, ...
        time_rife(i)*1000, time_mrife(i)*1000, time_irife(i)*1000, time_iirife(i)*1000);
end

% 分析计算复杂度
fprintf('\n=== 性能分析 ===\n');

% 找到最快和最慢的算法
times_all = [time_fft; time_czt; time_improved_czt; time_rife; time_mrife; time_irife; time_iirife];
algorithm_names = {'FFT-Peak', 'CZT', '改进CZT', 'RIFE', 'MRIFE', 'IRIFE', 'IIRIFE'};

% 对于最大采样点数的情况
max_idx = num_n;
[min_time, fastest_idx] = min(times_all(:, max_idx));
[max_time, slowest_idx] = max(times_all(:, max_idx));

fprintf('N = %d 时的性能对比:\n', N_values(max_idx));
fprintf('最快算法: %s (%.3f ms)\n', algorithm_names{fastest_idx}, min_time*1000);
fprintf('最慢算法: %s (%.3f ms)\n', algorithm_names{slowest_idx}, max_time*1000);
fprintf('速度差异: %.1f倍\n', max_time/min_time);

% 分析计算复杂度增长趋势
fprintf('\n计算复杂度增长分析 (N从256到16384):\n');
for alg = 1:length(algorithm_names)
    ratio = times_all(alg, end) / times_all(alg, 1);
    growth_factor = log2(ratio) / log2(N_values(end)/N_values(1));
    fprintf('%s: %.1f倍增长 (复杂度 ~O(N^%.2f))\n', algorithm_names{alg}, ratio, growth_factor);
end

% 相对性能比较
fprintf('\n相对性能比较 (以FFT为基准=1.0):\n');
for i = 1:num_n
    fprintf('N=%5d: ', N_values(i));
    for alg = 1:length(algorithm_names)
        relative_time = times_all(alg, i) / time_fft(i);
        fprintf('%s=%.1f ', algorithm_names{alg}(1:3), relative_time);
    end
    fprintf('\n');
end
