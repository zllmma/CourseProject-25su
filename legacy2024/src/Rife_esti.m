function [esti_freq, P1] = Rife_esti(x, t)
    L = length(t); %计算信号长度
    Fs = (L - 1) / (t(L) - t(1)); %计算采样频率

    % 计算DFT
    [~, P1, ~] = DFT(x, t);

    % 找到最大峰值的位置
    [~, idx] = max(P1);

    % 修正索引，如果峰值在第一个或最后一个点
    if idx == 1
        idx = 2;
    elseif idx == L / 2 + 1
        idx = L / 2;
    end

    % Rife算法
    if P1(idx + 1) <= P1(idx - 1)
        r = -1;
    else
        r = 1;
    end

    delta = P1(idx + r) / (P1(idx) + P1(idx + r));
    correction = r * delta * Fs / L;

    % 估计的频率
    esti_freq = (idx - 1) * Fs / L + correction;
end
