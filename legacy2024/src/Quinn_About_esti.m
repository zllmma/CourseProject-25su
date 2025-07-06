function [esti_freq, P1] = Quinn_About_esti(x, t)
    L = length(t); %计算信号长度
    Fs = (L - 1) / (t(L) - t(1)); %计算采样频率
    [X, P1, ~] = DFT(x, t);
    [~, idx] = max(P1);
    beta1 = real(X(idx - 1) / X(idx));
    beta2 = real(X(idx + 1) / X(idx));
    delta1 = beta1 / (1 - beta1);
    delta2 = beta2 / (beta2 - 1);

    if delta1 > 0 && delta2 > 0
        delta = delta2;
    else
        delta = delta1;
    end

    i = sqrt(-1);
    X05 = 0;

    for l = 1:L
        X05 = X05 + x(l) * exp(-i * (l - 1) * 2 * pi * (idx - 1 + 0.5) / L);
    end

    X_05 = 0;

    for l = 1:L
        X_05 = X_05 + x(l) * exp(-i * (l - 1) * 2 * pi * (idx - 0.5 - 1) / L);
    end

    X0 = 0;

    for l = 1:L
        X0 = X0 + x(l) * exp(-i * (l - 1) * 2 * pi * (idx \- 1) / L);
    end

    a = abs(X_05) / abs(X0);
    b = abs(X05) / abs(X0);

    if a > 1 || b > 1
        delta = 0.5 * real((X05 + X_05) / (X05 - X_05));
    end

    esti_freq = (idx - 1) * Fs / L + delta * Fs / L;
end
