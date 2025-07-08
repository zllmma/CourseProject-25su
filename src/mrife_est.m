function f_est = mrife_est(s, fs)
%   使用M-Rife算法精确估计信号频率
%   输入参数:
%       s  : 输入信号向量 (1 x N)
%       fs : 采样频率 (Hz)
%   输出参数:
%       f_est  : 估计的频率 (Hz)

% ===== 步骤1: 执行第一次Rife估计 =====
N = length(s);
[f1, ~, r1, delta_mag1] = rife_algorithm(s, fs);  % 调用Rife函数

% ===== 步骤2: 计算归一化频偏 =====
delta_f = r1 * delta_mag1 * (fs / N);  % 实际频偏
alpha = abs(delta_f) * N / fs;        % 归一化频偏 (0~0.5)

% ===== 步骤3: 判断是否需要进行频谱平移 =====
if alpha < 1/3
    % 如果在中心1/3区域，直接返回Rife估计结果
    f_est = f1;
else
    % ===== 步骤4: 计算频谱平移量 =====
    delta_shift = -sign(delta_f) * (fs / (3 * N));  % 论文公式(7)
    
    % ===== 步骤5: 生成频移信号 =====
    n = (0:N-1);
    shifted_signal = s .* exp(1j * 2 * pi * delta_shift * n / fs);  % 论文公式(8)
    
    % ===== 步骤6: 对频移信号执行第二次Rife估计 =====
    f2 = rife_algorithm(shifted_signal, fs);
    
    % ===== 步骤7: 计算最终频率估计 =====
    f_est = f2 - delta_shift;  % 论文公式(9)
end


