function [esti_freq, P1] = AIrife_esti(x, t)
    L = length(t); %计算信号长度
    Fs = (L - 1) / (t(L) - t(1)); %计算采样频率

    % 计算DFT
    [~, P1, ~] = DFT(x, t);

    % 找到最大峰值的位置
    [~, idx] = max(P1);

    fe = (idx - 1) * Fs / L;

    % 简单的频谱细化技术（线性插值）
    %alpha=simple_zfft(idx,L,P1);

    % 通过Zoom-FFT算法细化FFT %待完善
    [y, freq] = zoom_fft(x, fe, Fs);

    [~, q] = find(freq == fe);

    if abs(y(q + 1)) > abs(y(q - 1))
        a = 1;
    else
        a = -1;
    end

    if abs(y(q)) <= abs(y(q + a))
        delta_k = 0.5 - abs(y(q + 2 * a)) / (abs(y(q)) + abs(y(q + 2 * a)));
    else
        delta_k = abs(y(q - a)) / (abs(y(q - a)) + abs(y(q + a)));
    end

    i = sqrt(-1);
    X1 = 0;

    for l = 1:L
        X1 = X1 + x(l) * exp(-i * (l - 1) * 2 * pi * (idx - a * delta_k + a - 1) / L);
    end

    X2 = 0;

    for l = 1:L
        X2 = X2 + x(l) * exp(-i * (l - 1) * 2 * pi * (idx - a * delta_k - 1) / L);
    end

    alpha = a * abs(X1) / (abs(X2) + abs(X1)) - a * delta_k;

    esti_freq = fe + alpha * Fs / L;

end
