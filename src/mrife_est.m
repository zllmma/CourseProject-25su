function f_est = mrife_est(s, fs)
% M_RIFE_ALGORITHM 使用M-Rife算法精确估计信号频率
%   输入参数:
%       s  : 输入信号向量 (1 x N)
%       fs : 采样频率 (Hz)
%   输出参数:
%       f_est  : 估计的频率 (Hz)

N = length(s);  % 信号长度
% ===== 步骤1: 执行第一次Rife估计 =====
[f1, k1, ~, ~] = rife_algorithm(s, fs);

% ===== 步骤2: 计算归一化频偏 =====
bin_width = fs / N;  % FFT频率分辨率
f_m = (k1 - 1) * bin_width;  % 最大谱线对应的频率
delta_f = f1 - f_m;          % 实际频偏

% 计算归一化频偏 (以bin为单位)
delta_norm = abs(delta_f) / bin_width;

% ===== 步骤3: 判断是否需要进行频谱平移 =====
if delta_norm <= 1/3
    % 如果在中心1/3区域，直接返回Rife估计结果
    f_est = f1;
    return;
else
    % ===== 步骤4: 计算频谱平移量 =====
    % 根据公式(7)计算δ
    if delta_norm <= bin_width / (3 * bin_width)  % 原文条件：f_Rife - mf_s/n ≤ f_s/(3n)
        % 计算频谱幅度
        Y = fft(s, N);
        A_m = abs(Y(k1));
        
        % 确定相邻谱线位置
        if delta_f < 0
            k_adj = k1 - 1;
            if k_adj < 1, k_adj = N; end
        else
            k_adj = k1 + 1;
            if k_adj > N, k_adj = 1; end
        end
        A_adj = abs(Y(k_adj));
        
        % 计算δ
        delta_val = 0.5 - A_adj / (A_adj + A_m);
    else
        delta_val = 0;
    end
    
    % 确定平移方向
    r_sign = sign(delta_f);
    
    % 计算平移量 (根据公式(9)中的rδf_s/N)
    delta_shift = r_sign * delta_val * bin_width;
    
    % ===== 步骤5: 生成频移信号 =====
    n = (0:N-1);
    % 根据公式(8): x'(n) = x(n) * e^{jr^2π(σ/N)n}
    % 注意：原文中的σ对应这里的delta_shift，但公式(8)有笔误，应为j2π而不是jr^2π
    shifted_signal = s .* exp(1j * 2 * pi * delta_shift * n / fs);
    
    % ===== 步骤6: 对频移信号执行第二次Rife估计 =====
    [f2, ~, ~, ~] = rife_algorithm(shifted_signal, fs);
    
    % ===== 步骤7: 计算最终频率估计 =====
    % 根据公式(9): f_{M-Rife} = f_0' - rδf_s/N
    f_est = f2 - delta_shift;
end
end