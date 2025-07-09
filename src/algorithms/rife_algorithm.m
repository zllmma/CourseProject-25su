function [f0_hat, k0, r, delta_mag] = rife_algorithm(s, fs)
    %   使用Rife算法精确估计信号频率
    %   输入参数:
    %       s  : 输入信号向量 (1 x N)
    %       fs : 采样频率 (Hz)
    %   输出参数:
    %       f0_hat  : 估计的频率 (Hz)
    %       k0     : 最大谱线索引
    %       r      : 修正方向 (1 或 -1)
    %       delta_mag : 修正因子大小

    % 步骤1: 计算N点DFT
    N = length(s);
    S = fft(s, N);

    % 步骤2: 计算幅度谱
    mag = abs(S);

    % 步骤3: 找到最大谱线索引k0 (MATLAB索引从1开始)
    [~, k0] = max(mag);
    k0 = k0 - 1; % 转换为从0开始的索引 (0到N-1)

    % 步骤4: 处理边界情况 (循环频谱)
    k0_left = k0 - 1;

    if k0_left < 0
        k0_left = N - 1; % 循环到频谱末端
    end

    k0_right = k0 + 1;

    if k0_right >= N
        k0_right = 0; % 循环到频谱起始
    end

    % 步骤5: 确定修正方向r
    if mag(k0_right + 1) >= mag(k0_left + 1) % +1因MATLAB索引
        r = 1;
        k_second = k0_right; % 次大谱线索引
    else
        r = -1;
        k_second = k0_left; % 次大谱线索引
    end

    % 步骤6: 计算修正因子大小
    S_k0 = mag(k0 + 1); % 最大谱线幅值
    S_k_second = mag(k_second + 1); % 次大谱线幅值
    delta_mag = S_k_second / (S_k_second + S_k0);

    % 步骤7: 计算修正后的频率估计
    f0_hat = (fs / N) * (k0 + r * delta_mag);
end
