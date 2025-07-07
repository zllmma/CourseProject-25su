function f_c = improved_czt_est(s, fs, q, M)
% improved_czt_freq_est: 使用改进的 CZT 方法估计频率，该方法可以校正栅栏效应
%
% 输入:
%   s   - 输入信号向量
%   fs  - 采样频率 (Hz)
%   q   - 用于控制细化区间半宽度的整数
%   M   - 频率细化倍数 (CZT频谱中的点数)
%
% 输出:
%   f_c - 高精度估计频率 (Hz)
%
% 该函数在 CZT 方法的基础上，利用精细频谱中峰值及其相邻谱线的幅度来估计和校正频率偏差 (delta)

% 获取采样点数 N
N = length(s);

% --- 步骤 1: 使用 FFT 进行粗略频率估计 ---
S = fft(s);
[~, m0_matlab] = max(abs(S(1:floor(N/2))));
m0 = m0_matlab - 1;
delta_f0 = fs / N;

% --- 步骤 2: 使用 CZT 生成精细频谱 ---
theta0 = 2 * pi * (m0 - q) / N;
phi0 = 2 * pi * (2 * q) / (M * N);
A = exp(1j * theta0);
W = exp(-1j * phi0);
Sczt = czt(s, M, W, A);

% 搜索最大谱线及其索引 m1
[~, m1_matlab] = max(abs(Sczt));
m1 = m1_matlab - 1;

% --- 步骤 3: 频率校正 ---

% 处理边界情况：如果峰值位于 CZT 频谱的开始或结束处，则无法进行三点插值
% 在这种情况下，我们回退到标准的 CZT 估计
if m1 == 0 || m1 == M - 1
    % 回退到标准 CZT 估计
    delta_f1 = (2 * q * fs) / (M * N);
    f_start = (m0 - q) * delta_f0;
    f_c = f_start + m1 * delta_f1;
    return;
end

% 获取最大谱线及其直接相邻谱线的幅度
amp_m1_minus_1 = abs(Sczt(m1_matlab - 1));
amp_m1 = abs(Sczt(m1_matlab));
amp_m1_plus_1 = abs(Sczt(m1_matlab + 1));

% 计算幅度比 a1 和 a2
a1 = amp_m1_plus_1 / amp_m1;
a2 = amp_m1_minus_1 / amp_m1;

% 计算频率偏移量 delta (δ)
denominator_cos_term = cos(2 * q * pi / M);
delta = (a1 - a2) / (a1 + a2 - 2 * denominator_cos_term);

% --- 步骤 4: 计算最终的高精度频率 ---

% 真实频率 fc 是由 delta 校正后的 CZT 估计值
delta_f1 = (2 * q * fs) / (M * N);
f_start = (m0 - q) * delta_f0;

% 校正后的频率为 fc = f_start + (m1 + delta) * delta_f1
f_c = f_start + (m1 + delta) * delta_f1;

end