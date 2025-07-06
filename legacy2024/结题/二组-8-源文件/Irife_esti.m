function [esti_freq, P1, freq]  = Irife_esti(x,t)%esti_freq为估计出的频率，P1为DFT频谱，freq为细化后DFT频谱，x为待测信号，t为时间信号
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

    fe = (idx-1)*Fs/L; %直接估计得出的频率，考虑matlab向量起点为1，idx要-1
    
    % 简单的频谱细化技术（线性插值）（较简单） 
    %alpha=simple_zfft(idx,L,P1); 

    % 通过Zoom-FFT算法细化FFT(已预先设置好参数，是的能取到idx±0.5的频点)
    [y,freq] = zoom_fft(x,fe,Fs);

    %以下是I-Rife算法
    [~,q] = find(freq == fe); %确定细化后频谱原最大谱线的位置

    if abs(y(q+1)) > abs(y(q-1)) %计算频率修正方向
        a = 1;
    else
        a = -1;
    end

    delta_k = 0.5 - abs(y(q+2*a))/(abs(y(q))+abs(y(q+2*a))); %计算频率修正因子

    %单独计算频移后的两谱线
    i = sqrt(-1);%虚数i，防止冲突
    X1 = 0; %初始值
    for l = 1:L
       X1 = X1 + x(l)*exp(-i*(l-1)*2*pi*(idx-a*delta_k+a-1)/L); %考虑matlab向量起点为1，idx要-1
    end
    X2 = 0; %初始值
    for l = 1:L
       X2 = X2 + x(l)*exp(-i*(l-1)*2*pi*(idx-a*delta_k-1)/L); %考虑matlab向量起点为1，idx要-1
    end
    
    alpha = a*abs(X1)/(abs(X2)+abs(X1))-a*delta_k; %频率修正项

    esti_freq = fe + alpha*Fs/L;
    
end