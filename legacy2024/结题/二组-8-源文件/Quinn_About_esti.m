function [esti_freq, P1] = Quinn_About_esti(x,t)%esti_freq为估计出的频率，P1为DFT频谱，x为待测信号，t为时间信号
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

    %Quinn算法
    beta1 = real(X(idx-1)/X(idx)); %加入相位信息，取实部
    beta2 = real(X(idx+1)/X(idx));
    delta1 = beta1/(1-beta1);
    delta2 = beta2/(beta2-1);
    if delta1>0 && delta2>0 %beta为频率修正项
        delta = delta2;
    else
        delta = delta1;
    end
    %改进部分
    %单独计算频移后的三谱线
    i = sqrt(-1); %虚数i，防止冲突
    X05 = 0;
    for l = 1:L
       X05 = X05 + x(l)*exp(-i*(l-1)*2*pi*(idx-1+0.5)/L); %考虑matlab向量起点为1，idx要-1
    end
    X_05 = 0;
    for l = 1:L
       X_05 = X_05 + x(l)*exp(-i*(l-1)*2*pi*(idx-0.5-1)/L); %考虑matlab向量起点为1，idx要-1
    end
    X0 = 0;
    for l = 1:L
       X0 = X0 + x(l)*exp(-i*(l-1)*2*pi*(idx\-1)/L); %考虑matlab向量起点为1，idx要-1
    end
    a = abs(X_05)/abs(X0);
    b = abs(X05)/abs(X0);
    if a>1 || b>1 %判断频率修正项的取值
        delta = 0.5*real((X05+X_05)/(X05-X_05)); 
    end
    esti_freq = (idx-1)*Fs/L+delta*Fs/L; %考虑matlab向量起点为1，idx要-1
end