function f0_hat = irife_est(s, fs)
%   使用改进的I-Rife算法精确估计信号频率
%   输入参数:
%       s  : 输入信号向量 (1 x N)
%       fs : 采样频率 (Hz)
%   输出参数:
%       f0_hat : 估计的频率 (Hz)

% 步骤1: 计算N点DFT并找到最大谱线位置
N = length(s);
S = fft(s, N);
mag = abs(S);
[~, k0_idx] = max(mag);
k0 = k0_idx - 1; % 转换为从0开始的索引 (0到N-1)

% 处理频谱循环特性
k0_plus_05 = mod(k0 + 0.5, N);
k0_minus_05 = mod(k0 - 0.5, N);

% 步骤2: 计算k0±0.5处的频谱幅值 (使用精确的DTFT计算)
n = (0:N-1);
S_k0_plus_05 = sum(s .* exp(-1j * 2 * pi * k0_plus_05 * n / N));
S_k0_minus_05 = sum(s .* exp(-1j * 2 * pi * k0_minus_05 * n / N));

% 步骤3: 确定修正方向r (根据文献[15-16])
if abs(S_k0_plus_05) >= abs(S_k0_minus_05)
    r = 1;
else
    r = -1;
end

% 步骤4: 计算初始修正因子δ和频移因子Δk
% 计算相邻谱线位置 (考虑循环边界)
k_adj = mod(k0 + r, N);
k_adj_idx = k_adj + 1; % MATLAB索引

% 获取幅值
S_k0 = abs(S(k0_idx));
S_adj = abs(S(k_adj_idx));

% 关键修正1: 正确的δ计算公式
delta = S_adj / (S_adj + S_k0); % |δ|
shift_k = 0.5 - delta;          % Δk = 0.5 - δ

% 步骤5: 执行频移操作
shift_amount = r * shift_k;
s_shifted = s .* exp(-1j * 2 * pi * shift_amount * n / N);

% 步骤6: 计算频移后信号在目标位置的频谱幅值
% 目标位置 (考虑循环边界)
k1 = mod(k0 - r * shift_k, N);
k2 = mod(k1 + r, N);

% 计算DTFT幅值
S_k1 = sum(s_shifted .* exp(-1j * 2 * pi * k1 * n / N));
S_k2 = sum(s_shifted .* exp(-1j * 2 * pi * k2 * n / N));

% 关键修正2: 确保分母不为零
denom = abs(S_k1) + abs(S_k2);
if denom < eps
    delta_new = 0.5; % 默认值
else
    delta_new = abs(S_k2) / denom; % |δ_new|
end

% 步骤7: 计算最终频率估计
% 关键修正3: 使用原始k0而不是频移后位置
k_est = k0 - r * shift_k + r * delta_new;
f0_hat = (fs / N) * k_est;

% 关键修正4: 频率补偿
f0_hat = f0_hat + (r * shift_k * fs / N);

% 确保频率在[0, fs)范围内
f0_hat = mod(f0_hat, fs);
end