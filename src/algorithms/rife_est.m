function f0_hat = rife_est(s, fs)
    %   使用Rife算法精确估计信号频率
    %   输入参数:
    %       s  : 输入信号向量 (1 x N)
    %       fs : 采样频率 (Hz)
    %   输出参数:
    %       f0_hat  : 估计的频率 (Hz)

    [f0_hat, ~, ~, ~] = rife_algorithm(s, fs);
end
