function [esti_freq, P1] = Quinn_esti(x,t)%esti_freq为估计出的频率，P1为DFT频谱，x为待测信号，t为时间信号
    L = length(t); %计算信号长度
    Fs = (L-1)/(t(L)-t(1)); %计算采样频率
    [X,P1,~] = DFT(x,t); %做FFT，得出DFT频谱
    [~,idx] = max(P1); %搜索找到最大谱线位置

    % 修正索引，如果峰值在第一个或最后一个点  
    if idx == 1  
        idx = 2;  
    elseif idx == L/2+1  
        idx = L/2;  
    end  
    
    %以下是Quinn算法
    beta1 = real(X(idx-1)/X(idx)); %加入相位信息，取实部
    beta2 = real(X(idx+1)/X(idx));
    delta1 = beta1/(1-beta1);
    delta2 = beta2/(beta2-1);
    if delta1>0 && delta2>0 %beta为频率修正项
        delta = delta2;
    else
        delta = delta1;
    end
    esti_freq = (idx-1)*Fs/L+delta*Fs/L; %考虑matlab向量起点为1，idx要-1
end