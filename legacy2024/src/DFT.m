function [X, P1, P2] = DFT(x, t)
    % 计算DFT
    L = length(t); %计算信号长度
    X = fft(x);
    P2 = abs(X / L);
    P1 = P2(1:L / 2 + 1);
    P1(2:end - 1) = 2 * P1(2:end - 1);
end
