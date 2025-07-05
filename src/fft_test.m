fs = 1000;         % 采样频率
N = 1024;           % 采样点数
t = (0:N-1)/fs;    % 时间向量
f = 50;            % 频率
x = sin(2 * pi * f * t);
Y = fft(x);
F = fs / N;        % 频率分辨率
f_fft = F*(0: N - 1);
plot(f_fft, abs(Y) / fs);