function [esti_freq, P1] = Quinn_esti(x, t)
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

    esti_freq = (idx - 1) * Fs / L + delta * Fs / L;
end
