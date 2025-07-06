function [esti_freq, P1] = Dirc_esti(x, t)
    L = length(t); %计算信号长度
    Fs = (L - 1) / (t(L) - t(1)); %计算采样频率

    % 计算DFT
    [~, P1, ~] = DFT(x, t);

    % 找到最大峰值的位置
    [~, idx] = max(P1);
    esti_freq = (idx - 1) * Fs / L;
end
