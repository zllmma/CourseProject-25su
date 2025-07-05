% 示例：使用CZT进行频率估计
fs = 1000;         % 采样频率
N = 256;           % 采样点数
t = (0:N-1)/fs;    % 时间向量

f_true = 50.3;     % 真实信号频率 (故意设置一个不在FFT谱线上的频率)
x = sin(2*pi*f_true*t) + 0.1*randn(size(t)); % 含噪声的正弦信号

% 1. 粗略FFT分析
Y_fft = fft(x);
P2 = abs(Y_fft/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f_fft = fs*(0:(N/2))/N;

figure;
subplot(2,1,1);
plot(f_fft, P1);
title('FFT 频谱 (粗略)');
xlabel('频率 (Hz)');
ylabel('幅度');
grid on;

% 2. 确定感兴趣的频率范围 (例如，在FFT峰值附近)
% 假设我们从FFT看出峰值在50Hz左右
f_start = 49;      % CZT起始频率
f_end = 51;        % CZT终止频率
M_czt = 2000;      % CZT输出点数 (越多分辨率越高)

% 使用 czt 函数进行频率估计
% cz(x, M, W, A)
% x: 输入信号
% M: 输出点数
% W: Z平面上的螺旋线的比率 W = exp(-j*2*pi*(f_end - f_start) / (M * fs))
% A: Z平面上的螺旋线的起始点 A = exp(j*2*pi*f_start / fs)
% 这里我们关心的是频率轴，所以W和A的设置与频率范围相关

% 转换到 C. Z. Transform 参数
W = exp(-1i * 2 * pi * (f_end - f_start) / (M_czt * fs));
A = exp(1i * 2 * pi * f_start / fs);

Y_czt = czt(x, M_czt, W, A);

% 计算CZT对应的频率轴
f_czt = (f_start + (0:M_czt-1) * (f_end - f_start) / M_czt);

subplot(2,1,2);
plot(f_czt, abs(Y_czt));
title('CZT 频谱 (局部放大)');
xlabel('频率 (Hz)');
ylabel('幅度');
grid on;

% 寻找CZT频谱的峰值
[max_val, idx_max] = max(abs(Y_czt));
estimated_f_czt = f_czt(idx_max);

fprintf('FFT 估计的频率 (粗略): %.2f Hz\n', f_fft(find(P1==max(P1))));
fprintf('CZT 估计的频率 (精确): %.4f Hz\n', estimated_f_czt);
fprintf('真实频率: %.2f Hz\n', f_true);