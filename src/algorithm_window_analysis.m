% =========================================================================
%     不同算法在不同窗函数和高低信噪比下的RMSE和标准差对比分析
% =========================================================================
%
% 功能说明:
%   1. 对比七种频率估计算法在不同窗函数下的性能
%   2. 分析高信噪比(20dB)和低信噪比(0dB)条件下的差异
%   3. 评估RMSE和标准差两个关键指标
%   4. 生成详细的对比图表和分析报告
%
% 算法列表: FFT-Peak, CZT, 改进CZT, Rife, M-Rife, I-Rife, IIRife
% 窗函数列表: 无窗, Hann窗, Hamming窗, Blackman窗, Kaiser窗
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

% --- 1. 仿真参数设置 ---
fs = 200e6;                    % 采样频率 (200 MHz)
N = 1024;                      % 采样点数
t = (0:N-1) / fs;             % 时间向量
A = 1.0;                      % 信号幅度
f_center = 50e6;              % 中心频率 (50 MHz)
delta_f0 = fs / N;            % 频率分辨率

% 高低信噪比设置
SNR_high = 20;                % 高信噪比 (dB)
SNR_low = 0;                  % 低信噪比 (dB)

% 相对频偏范围
relative_offsets = -0.5:0.1:0.5;  % 减少点数以加快仿真
num_offsets = length(relative_offsets);

num_trials = 1000;            % 蒙特卡洛试验次数

% CZT 参数
q = 1;
M = 64;

% 窗函数配置 (移除矩形窗)
window_types = {'none', 'hann', 'hamming', 'blackman', 'kaiser'};
window_names = {'无窗', 'Hann窗', 'Hamming窗', 'Blackman窗', 'Kaiser窗'};
num_windows = length(window_types);

% 算法名称
algorithm_names = {'FFT-Peak', 'CZT', '改进CZT', 'Rife', 'M-Rife', 'I-Rife', 'IIRife'};
num_algorithms = length(algorithm_names);

% --- 2. 初始化结果存储 ---
% 维度：[窗函数, 信噪比(1=低,2=高), 频偏, 算法]
rmse_results = zeros(num_windows, 2, num_offsets, num_algorithms);
std_results = zeros(num_windows, 2, num_offsets, num_algorithms);

fprintf('=== 开始算法在不同窗函数下的性能对比分析 ===\n');
fprintf('窗函数数量: %d (无窗, Hann, Hamming, Blackman, Kaiser), 算法数量: %d\n', num_windows, num_algorithms);
fprintf('高SNR: %d dB, 低SNR: %d dB\n', SNR_high, SNR_low);
fprintf('频偏点数: %d, 试验次数: %d\n', num_offsets, num_trials);

% --- 3. 主分析循环 ---
for w_idx = 1:num_windows
    window_type = window_types{w_idx};
    window_name = window_names{w_idx};
    
    fprintf('\n--- 处理窗函数: %s ---\n', window_name);
    % 生成窗函数 (不包含矩形窗)
    switch window_type
        case 'none'
            w = ones(1, N);
        case 'hann'
            w = hann(N)';
        case 'hamming'
            w = hamming(N)';
        case 'blackman'
            w = blackman(N)';
        case 'kaiser'
            beta = 8.6;
            w = kaiser(N, beta)';
    end
    
    % 窗函数归一化
    if ~strcmp(window_type, 'none')
        w = w / sqrt(mean(w.^2));
    end
    
    % 对高低两种信噪比进行分析
    for snr_idx = 1:2
        if snr_idx == 1
            SNR_dB = SNR_low;
            snr_label = '低SNR';
        else
            SNR_dB = SNR_high;
            snr_label = '高SNR';
        end
        
        fprintf('  %s (%d dB) 分析中...\n', snr_label, SNR_dB);
        
        % 计算噪声参数
        snr_linear = 10^(SNR_dB/10);
        
        % 对每个频偏进行分析
        for off_idx = 1:num_offsets
            current_offset = relative_offsets(off_idx);
            f_true = f_center + current_offset * delta_f0;
            
            fprintf('    频偏 %.1f 处理中...\n', current_offset);
            
            % 存储所有算法在当前条件下的估计误差
            errors_all = zeros(num_algorithms, num_trials);
            
            % 蒙特卡洛仿真
            for trial = 1:num_trials
                % 生成纯净信号
                phi = 2 * pi * rand();
                s_clean = A * exp(1j * (2*pi*f_true*t + phi));
                
                % 应用窗函数
                s_windowed = s_clean .* w;
                
                % 计算实际噪声功率
                windowed_signal_power = mean(abs(s_windowed).^2);
                actual_snr_linear = 10^(SNR_dB/10);
                actual_noise_power = windowed_signal_power / actual_snr_linear;
                actual_noise_std = sqrt(actual_noise_power / 2);
                
                % 生成噪声
                noise = (randn(1,N) + 1j*randn(1,N)) * actual_noise_std;
                noise_windowed = noise .* w;
                s_noisy = s_windowed + noise_windowed;
                
                % 执行所有算法
                try
                    f_estimates = [
                        fft_est(s_noisy, fs),...
                        czt_est(s_noisy, fs, q, M),...
                        improved_czt_est(s_noisy, fs, q, M),...
                        rife_est(s_noisy, fs),...
                        mrife_est(s_noisy, fs),...
                        irife_est(s_noisy, fs),...
                        iirife_est(s_noisy, fs)
                    ];
                    
                    % 计算误差
                    errors_all(:, trial) = f_estimates' - f_true;
                catch
                    errors_all(:, trial) = NaN;
                end
            end
            
            % 计算RMSE和标准差
            for alg_idx = 1:num_algorithms
                % 直接计算误差的RMSE和标准差，忽略NaN值
                valid_errors = errors_all(alg_idx, ~isnan(errors_all(alg_idx, :)));
                rmse_results(w_idx, snr_idx, off_idx, alg_idx) = sqrt(mean(valid_errors.^2));
                std_results(w_idx, snr_idx, off_idx, alg_idx) = std(valid_errors);
            end

        end
    end
end

fprintf('\n仿真完成！开始生成分析结果...\n');

% --- 4. 数据分析 ---

% 计算平均性能指标（跨频偏平均）
avg_rmse = squeeze(mean(rmse_results, 3, 'omitnan'));  % [窗函数, 信噪比, 算法]
avg_std = squeeze(mean(std_results, 3, 'omitnan'));

% --- 5. 可视化结果 ---

% 图1: 算法在不同窗函数下的RMSE对比（低SNR）
figure;

subplot(2, 3, 1);
rmse_low_snr = squeeze(avg_rmse(:, 1, :))';  % [算法, 窗函数]
bar(rmse_low_snr);
set(gca, 'XTickLabel', algorithm_names, 'XTickLabelRotation', 45);
ylabel('平均RMSE (Hz)');
title(sprintf('低SNR (%d dB) - 算法RMSE对比', SNR_low));
legend(window_names, 'Location', 'best');
grid on;
set(gca, 'YScale', 'log');

% 图2: 算法在不同窗函数下的RMSE对比（高SNR）
subplot(2, 3, 2);
rmse_high_snr = squeeze(avg_rmse(:, 2, :))';  % [算法, 窗函数]
bar(rmse_high_snr);
set(gca, 'XTickLabel', algorithm_names, 'XTickLabelRotation', 45);
ylabel('平均RMSE (Hz)');
title(sprintf('高SNR (%d dB) - 算法RMSE对比', SNR_high));
legend(window_names, 'Location', 'best');
grid on;
set(gca, 'YScale', 'log');

% 图3: 算法在不同窗函数下的标准差对比（低SNR）
subplot(2, 3, 3);
std_low_snr = squeeze(avg_std(:, 1, :))';  % [算法, 窗函数]
bar(std_low_snr);
set(gca, 'XTickLabel', algorithm_names, 'XTickLabelRotation', 45);
ylabel('平均标准差 (Hz)');
title(sprintf('低SNR (%d dB) - 算法标准差对比', SNR_low));
legend(window_names, 'Location', 'best');
grid on;
set(gca, 'YScale', 'log');

% 图4: 算法在不同窗函数下的标准差对比（高SNR）
subplot(2, 3, 4);
std_high_snr = squeeze(avg_std(:, 2, :))';  % [算法, 窗函数]
bar(std_high_snr);
set(gca, 'XTickLabel', algorithm_names, 'XTickLabelRotation', 45);
ylabel('平均标准差 (Hz)');
title(sprintf('高SNR (%d dB) - 算法标准差对比', SNR_high));
legend(window_names, 'Location', 'best');
grid on;
set(gca, 'YScale', 'log');

% 图5: SNR改善比率（RMSE）
subplot(2, 3, 5);
rmse_improvement = squeeze(avg_rmse(:, 1, :)) ./ squeeze(avg_rmse(:, 2, :));  % 低SNR/高SNR
rmse_improvement = rmse_improvement';  % [算法, 窗函数]
bar(rmse_improvement);
set(gca, 'XTickLabel', algorithm_names, 'XTickLabelRotation', 45);
ylabel('RMSE改善比率 (低SNR/高SNR)');
title('算法RMSE随SNR的改善程度');
legend(window_names, 'Location', 'best');
grid on;

% 图6: SNR改善比率（标准差）
subplot(2, 3, 6);
std_improvement = squeeze(avg_std(:, 1, :)) ./ squeeze(avg_std(:, 2, :));  % 低SNR/高SNR
std_improvement = std_improvement';  % [算法, 窗函数]
bar(std_improvement);
set(gca, 'XTickLabel', algorithm_names, 'XTickLabelRotation', 45);
ylabel('标准差改善比率 (低SNR/高SNR)');
title('算法标准差随SNR的改善程度');
legend(window_names, 'Location', 'best');
grid on;

sgtitle('不同算法在各种窗函数和高低信噪比条件下的性能对比', 'FontSize', 14);

% --- 6. 详细的频偏特性分析图 ---
figure;

% 选择几个代表性算法进行详细分析
representative_algs = [1, 3, 6];  % FFT-Peak, 改进CZT, I-Rife
alg_titles = {'FFT-Peak算法', '改进CZT算法', 'I-Rife算法'};

for i = 1:3
    alg_idx = representative_algs(i);
    
    % RMSE随频偏变化（低SNR）
    subplot(3, 4, (i-1)*4 + 1);
    for w_idx = 1:num_windows
        plot(relative_offsets, squeeze(rmse_results(w_idx, 1, :, alg_idx)), ...
             '-o', 'LineWidth', 1.5, 'DisplayName', window_names{w_idx});
        hold on;
    end
    hold off;
    grid on;
    set(gca, 'YScale', 'log');
    xlabel('相对频偏');
    ylabel('RMSE (Hz)');
    title([alg_titles{i} ' - 低SNR RMSE']);
    if i == 1, legend('show', 'Location', 'best'); end
    
    % RMSE随频偏变化（高SNR）
    subplot(3, 4, (i-1)*4 + 2);
    for w_idx = 1:num_windows
        plot(relative_offsets, squeeze(rmse_results(w_idx, 2, :, alg_idx)), ...
             '-o', 'LineWidth', 1.5, 'DisplayName', window_names{w_idx});
        hold on;
    end
    hold off;
    grid on;
    set(gca, 'YScale', 'log');
    xlabel('相对频偏');
    ylabel('RMSE (Hz)');
    title([alg_titles{i} ' - 高SNR RMSE']);
    
    % 标准差随频偏变化（低SNR）
    subplot(3, 4, (i-1)*4 + 3);
    for w_idx = 1:num_windows
        plot(relative_offsets, squeeze(std_results(w_idx, 1, :, alg_idx)), ...
             '-o', 'LineWidth', 1.5, 'DisplayName', window_names{w_idx});
        hold on;
    end
    hold off;
    grid on;
    set(gca, 'YScale', 'log');
    xlabel('相对频偏');
    ylabel('标准差 (Hz)');
    title([alg_titles{i} ' - 低SNR 标准差']);
    
    % 标准差随频偏变化（高SNR）
    subplot(3, 4, (i-1)*4 + 4);
    for w_idx = 1:num_windows
        plot(relative_offsets, squeeze(std_results(w_idx, 2, :, alg_idx)), ...
             '-o', 'LineWidth', 1.5, 'DisplayName', window_names{w_idx});
        hold on;
    end
    hold off;
    grid on;
    set(gca, 'YScale', 'log');
    xlabel('相对频偏');
    ylabel('标准差 (Hz)');
    title([alg_titles{i} ' - 高SNR 标准差']);
end

sgtitle('代表性算法在不同窗函数下的频偏特性详细分析', 'FontSize', 14);

% --- 7. 生成分析报告 ---
fprintf('\n=== 算法性能分析报告 ===\n');

% 找出每种条件下的最佳算法
fprintf('\n1. 各种窗函数下的最佳算法（基于平均RMSE）:\n\n');

fprintf('%-12s | %-15s | %-15s\n', '窗函数', '低SNR最佳算法', '高SNR最佳算法');
fprintf(repmat('-', 1, 50));
fprintf('\n');

for w_idx = 1:num_windows
    [~, best_low] = min(avg_rmse(w_idx, 1, :));
    [~, best_high] = min(avg_rmse(w_idx, 2, :));
    fprintf('%-12s | %-15s | %-15s\n', window_names{w_idx}, ...
            algorithm_names{best_low}, algorithm_names{best_high});
end

% 算法排名分析
fprintf('\n2. 算法综合性能排名（所有窗函数平均）:\n\n');

% 计算每个算法的综合得分
overall_scores_low = squeeze(mean(avg_rmse(:, 1, :), 1));
overall_scores_high = squeeze(mean(avg_rmse(:, 2, :), 1));

[~, ranking_low] = sort(overall_scores_low);
[~, ranking_high] = sort(overall_scores_high);

fprintf('低SNR排名:\n');
for i = 1:num_algorithms
    alg_idx = ranking_low(i);
    fprintf('%d. %-12s (RMSE: %.2e Hz)\n', i, algorithm_names{alg_idx}, overall_scores_low(alg_idx));
end

fprintf('\n高SNR排名:\n');
for i = 1:num_algorithms
    alg_idx = ranking_high(i);
    fprintf('%d. %-12s (RMSE: %.2e Hz)\n', i, algorithm_names{alg_idx}, overall_scores_high(alg_idx));
end

% 窗函数敏感性分析
fprintf('\n3. 算法对窗函数的敏感性分析:\n\n');

for alg_idx = 1:num_algorithms
    % 计算该算法在不同窗函数下的性能变异系数
    cv_low = std(avg_rmse(:, 1, alg_idx)) / mean(avg_rmse(:, 1, alg_idx));
    cv_high = std(avg_rmse(:, 2, alg_idx)) / mean(avg_rmse(:, 2, alg_idx));
    
    fprintf('%-12s: 低SNR变异系数=%.3f, 高SNR变异系数=%.3f\n', ...
            algorithm_names{alg_idx}, cv_low, cv_high);
end

% 关键发现总结
fprintf('\n=== 关键发现 ===\n');

% 找出最不敏感的算法
cv_scores = zeros(num_algorithms, 2);
for alg_idx = 1:num_algorithms
    cv_scores(alg_idx, 1) = std(avg_rmse(:, 1, alg_idx)) / mean(avg_rmse(:, 1, alg_idx));
    cv_scores(alg_idx, 2) = std(avg_rmse(:, 2, alg_idx)) / mean(avg_rmse(:, 2, alg_idx));
end

[~, most_stable_low] = min(cv_scores(:, 1));
[~, most_stable_high] = min(cv_scores(:, 2));

fprintf('1. 低SNR下最稳定算法: %s (变异系数: %.3f)\n', ...
        algorithm_names{most_stable_low}, cv_scores(most_stable_low, 1));
fprintf('2. 高SNR下最稳定算法: %s (变异系数: %.3f)\n', ...
        algorithm_names{most_stable_high}, cv_scores(most_stable_high, 2));

fprintf('3. 最佳低SNR算法: %s\n', algorithm_names{ranking_low(1)});
fprintf('4. 最佳高SNR算法: %s\n', algorithm_names{ranking_high(1)});


fprintf('分析完成\n');
