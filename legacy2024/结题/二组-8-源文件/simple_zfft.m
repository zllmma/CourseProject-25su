function alpha=simple_zfft(idx,L,P)%简单插值细化FFT，用于改进Rife算法
    if idx > 1 && idx <= L/2
        leftidx = idx - 1;  
        rightidx = idx + 1;   
        % 确保索引在有效范围内  
        if leftidx < 1  
            leftidx = idx;  
        end  
        if rightidx > L/2  
            rightidx = idx;  
        end  

        alpha = (P(rightidx)-P(leftidx))/(2*P(idx)-(P(leftidx)+P(rightidx))/2); % 计算插值系数 
    else
        alpha = 0;
    end
end