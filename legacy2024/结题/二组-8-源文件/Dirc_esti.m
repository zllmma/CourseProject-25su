function [esti_freq, P1]= Dirc_esti(x,t)%esti_freq为估计出的频率，P1为DFT频谱，x为待测信号，t为时间信号
    L = length(t); %计算信号长度
    Fs = (L-1)/(t(L)-t(1)); %计算采样频率
    [~,P1,~]  = DFT(x,t); %做FFT，得出DFT频谱
    [~, idx] = max(P1); %搜索找到最大谱线位置
    
    % 修正索引，如果峰值在第一个或最后一个点  
    if idx == 1  
        idx = 2;  
    elseif idx == L/2+1  
        idx = L/2;  
    end  

    esti_freq = (idx-1)*Fs/L; %直接估计得出的频率，考虑matlab向量起点为1，idx要-1
end