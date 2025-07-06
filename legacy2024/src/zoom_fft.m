%%待完善
function [y, freq] = zoom_fft(x, fe, Fs) %细化FFT，用于改进Rife算法
    %x被测信号，fe是细化区间的中心频率，m是细化倍数
    m = 16;
    L = length(x); % 计算读入数据长度
    L_fft = 2 * L / m;
    fi = fe - Fs / m / 2; % 计算细化截止频率下限
    fa = fi + Fs / m; % 计算细化截止频率上限
    La = round(0.5 * L / m + 1); % 确定低通滤波器截止频率对应的谱线条数
    % 频移
    idx = 0:L - 1; % 序列索引号
    b = idx * pi * (fi + fa) / Fs; % 设置单位旋转因子
    y = x .* exp(-1i * b); % 进行频移
    X = fft(y, L); % FFT
    % 低通滤波和下采样
    a(1:La) = X(1:La); % 取正频率部分的低频成分
    a(L - La + 2:L) = X(L - La + 2:L); % 取负频率部分的低频成分
    X = ifft(a, L);
    c = X(1:m:L); % 下采样
    % 求细化频谱
    y = fft(c, L_fft) * 2 / L_fft; % 再一次FFT
    y = fftshift(y); % 重新排列
    freq = fi + (0:L_fft - 1) * Fs / m / L_fft; % 频率设置
end
