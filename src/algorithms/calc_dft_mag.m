function mag = calc_dft_mag(s, k_bin, N)
    % 计算任意分数索引k_bin处的DFT幅值
    %   s: 时域信号 (列向量)
    %   k_bin: 分数索引 (实数, 单位: bin)
    %   N: 信号长度
    k_bin = mod(k_bin, N); % 处理超出[0,N)范围的情况
    n = (0:N - 1);
    complex_exponent = -1j * 2 * pi * k_bin * n / N;
    dft_val = sum(s .* exp(complex_exponent));
    mag = abs(dft_val);
end
