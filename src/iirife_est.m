function f0_hat = iirife_est(s, fs)
%   基于插值修正的I-Rife频率估计算法
%   输入参数:
%       s  : 输入信号向量 (1 x N)
%       fs : 采样频率 (Hz)
%   输出参数:
%       f0_hat : 估计的频率 (Hz)

% 步骤1: 计算N点FFT并找到最大谱线索引
N = length(s);
S = fft(s, N);
mag = abs(S);
[~, idx] = max(mag);
k0 = idx - 1; % 转换为0-based索引

% 步骤2: 计算k0-0.5和k0+0.5处的幅值
k_minus_half = mod(k0 - 0.5, N);
k_plus_half = mod(k0 + 0.5, N);
mag_minus_half = calc_dft_mag(s, k_minus_half, N);
mag_plus_half = calc_dft_mag(s, k_plus_half, N);

% 步骤3: 确定修正方向r
if mag_plus_half > mag_minus_half
    r = 1;
else
    r = -1;
end

% 步骤4: 在修正方向插值k0+0.25r
k_interp = mod(k0 + 0.25 * r, N);
mag_interp = calc_dft_mag(s, k_interp, N);

% 步骤5: 获取三个关键点的幅值
mag_k0 = mag(idx); % 最大谱线幅值 (直接从FFT结果获取)

% 根据r选择k0+0.5r点的幅值
if r == 1
    mag_plus_half_r = mag_plus_half;
else
    mag_plus_half_r = mag_minus_half;
end

% 步骤6: 确定真实频率所在的小区域
if (mag_k0 >= mag_interp) && (mag_interp > mag_plus_half_r)
    region = 1; % (k0, k0+0.125r)
elseif (mag_plus_half_r >= mag_interp) && (mag_interp > mag_k0)
    region = 2; % (k0+0.375r, k0+0.5r)
else
    region = 3; % (k0+0.125r, k0+0.375r)
end

% 步骤7: 根据区域计算频移因子delta_k
switch region
    case 1 % (k0, k0+0.125r)
        % 使用|S(k0-0.5r)|和|S(k0+0.5r)|
        if r == 1
            mag_left = mag_minus_half;
            mag_right = mag_plus_half;
        else
            mag_left = mag_plus_half;
            mag_right = mag_minus_half;
        end
        delta_k = 0.5 * mag_left / (mag_left + mag_right);
        
    case 2 % (k0+0.375r, k0+0.5r)
        % 计算|S(k0+r)|
        k_plus_r = mod(k0 + r, N);
        mag_plus_r = calc_dft_mag(s, k_plus_r, N);
        delta_k = 0.5 * mag_plus_r / (mag_plus_r + mag_k0);
        
    case 3 % (k0+0.125r, k0+0.375r)
        % 使用|S(k0)|和|S(k0+0.5r)|
        delta_k = 0.5 * mag_k0 / (mag_k0 + mag_plus_half_r);
end

% 步骤8: 根据公式(10)计算最终频率估计
% 计算中间点 k_mid = k0 - r * delta_k
k_mid = k0 - r * delta_k;

% 计算两个DFT点的位置
k1 = mod(k_mid, N);            % k0 - r * delta_k
k2 = mod(k_mid + r, N);        % k0 - r * delta_k + r

% 计算这两个点的幅值
mag_k1 = calc_dft_mag(s, k1, N);
mag_k2 = calc_dft_mag(s, k2, N);

% 计算修正项
correction = r * mag_k2 / (mag_k1 + mag_k2);

% 最终频率估计
f0_hat = (fs / N) * (k_mid + correction);

end