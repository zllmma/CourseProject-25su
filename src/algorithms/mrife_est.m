function f_est = mrife_est(s, fs)
% MRIFE_EST 使用M-Rife算法精确估计信号频率
%   该函数实现了修正的Rife(MRife)算法，通过对信号进行频谱平移，
%   使其落在常规Rife算法性能最佳的区域，从而获得高精度且稳定的频率估计。
%
%   输入参数:
%       s  : 输入信号向量 (1 x N)，应为复数或实数信号
%       fs : 采样频率 (Hz)
%
%   输出参数:
%       f_est  : 估计的频率 (Hz)
%

% 步骤 1: 获取信号参数
N = length(s);
n = 0:(N-1); % 创建时间索引向量，从0到N-1

% 步骤 2: 对原始信号进行频谱平移
% 目标是平移 0.5 个频率分辨单元 (frequency bin)。
% 一个频率分辨单元的大小为 fs/N。
% 在时域上，这等同于将信号序列 x(n) 乘以一个复指数序列 exp(j*pi*n/N)。
% MATLAB 中虚数单位是 '1j' 或 '1i'。
shift_factor = exp(1j * pi * n / N);

% 产生经过频谱平移的新信号
s_shifted = s .* shift_factor; % 使用逐元素相乘

% 步骤 3: 对平移后的信号使用您提供的常规 Rife 算法
% 因为此时信号的频谱已经被移动到Rife算法性能最好的区域，
% 所以这里得到的估计结果非常精确。
% 我们只需要估计出的频率，所以用 ~ 忽略 rife_algorithm 的其他输出。
[f_shifted_est, ~, ~, ~] = rife_algorithm(s_shifted, fs);

% 步骤 4: 修正频率估计值
% 上一步得到的是“平移后”信号的频率，我们需要减去当初施加的平移量，
% 以此来还原出原始信号的真实频率。
freq_shift_amount = 0.5 * (fs / N); % 我们施加的平移量大小
f_est = f_shifted_est - freq_shift_amount;

end