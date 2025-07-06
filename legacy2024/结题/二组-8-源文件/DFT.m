function [X,P1,P2] = DFT(x,t)
    % 计算DFT  
    L = length(t); %计算信号长度
    X = fft(x); %fft
    P2 = abs(X/L); %计算谱线对应幅度大小
    P1 = P2(1:L/2+1); %由于对称，取前半
    P1(2:end-1) = 2*P1(2:end-1); %由于折叠，2倍
end