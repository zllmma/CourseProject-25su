function f_czt = czt_est(s, fs, q, M)
    % czt_freq_est: 使用 Chirp-Z 变换 (CZT) 方法估计频率
    %
    % 输入:
    %   s   - 输入信号向量
    %   fs  - 采样频率 (Hz)
    %   q   - 用于控制细化区间半宽度的整数
    %   M   - 频率细化倍数 (CZT频谱中的点数)
    %
    % 输出:
    %   f_czt - 估计的频率 (Hz)
    %
    % 该函数首先执行粗略的 FFT 以找到近似频率，然后使用 CZT 在该频率上进行放大以获得更精确的估计

    % 获取采样点数 N
    N = length(s);

    % --- 步骤 1: 使用 FFT 进行粗略频率估计 ---
    S = fft(s);
    [~, m0] = max(abs(S(1:floor(N / 2))));
    m0 = m0 - 1; % 最大谱线的索引
    delta_f0 = fs / N; % FFT 频率分辨率

    % --- 步骤 2: 使用 CZT 进行精细频率估计 ---

    % 定义 CZT 参数 A 和 W，用于 Z 平面单位圆上的螺旋轮廓
    % 轮廓从角度 theta0 开始，并扫过 phi0 * (M-1) 的角度

    % CZT 轮廓的起始角度
    theta0 = 2 * pi * (m0 - q) / N;

    % 相邻 CZT 采样点之间的角度步长
    phi0 = 2 * pi * (2 * q) / (M * N);

    % z_k = A * W^(-k)。对于单位圆缩放，设 A0 = 1, W0 = 1
    A = exp(1j * theta0);
    W = exp(-1j * phi0);

    % 对信号 s 执行 CZT 以获得精细频谱
    Sczt = czt(s, M, W, A);

    % 在 CZT 频谱中搜索峰值幅度
    [~, m1] = max(abs(Sczt));
    m1 = m1 - 1; % 精细频谱中最大峰值的索引

    % 计算 CZT 的频率分辨率
    delta_f1 = (2 * q * fs) / (M * N);

    % CZT 分析窗口的起始频率
    f_start = (m0 - q) * delta_f0;

    % 计算最终的、更精确的频率估计值
    f_czt = f_start + m1 * delta_f1;

end
