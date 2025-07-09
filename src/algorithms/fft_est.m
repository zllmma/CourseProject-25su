function f_est = fft_est(s, fs)
% fft_freq_est: 通过直接FFT谱峰检测法估计频率
%
% 输入:
%   s   - 输入信号向量
%   fs  - 采样频率 (Hz)
%
% 输出:
%   f_est - 估计的频率 (Hz)
%

% 获取采样点数 N
N = length(s);

% 执行 N 点离散傅里叶变换 (DFT)
S = fft(s);

% 在单边频谱中找到最大幅度的索引。
% 搜索范围限制在前半部分 (1 到 N/2) 以避免混叠
[~, m0] = max(abs(S(1:floor(N/2))));

% 公式中使用的索引 m0 是 MATLAB 索引减 1
m0 = m0 - 1;

% 计算 FFT 的频率分辨率
delta_f0 = fs / N;

% 通过将谱峰索引乘以分辨率来估计频率
f_est = m0 * delta_f0;

end