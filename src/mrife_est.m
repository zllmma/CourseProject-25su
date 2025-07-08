function f_est = mrife_est(s, fs)
%   智能选择Rife或M-Rife算法进行精确频率估计
%   该函数首先使用常规Rife算法进行探测，通过其修正因子delta_mag判断
%   信号频率是否落在Rife算法的“甜点区”。
%   如果落在甜点区，则直接采用Rife的结果；否则，启动更稳健的M-Rife算法。
%
%   输入参数:
%       s  : 输入信号向量 (1 x N)
%       fs : 采样频率 (Hz)
%
%   输出参数:
%       f_est  : 估计的频率 (Hz)
%

% 步骤 1: 进行一次“侦察性”的Rife算法计算
% 我们需要它的初步估计结果 f_rife，以及关键的判断指标 delta_mag。
[f_rife, ~, ~, delta_mag] = rife_algorithm(s, fs);


% 步骤 2: 设定判断阈值
% delta_mag 理论范围在 [0, 0.5] 之间 (考虑r的方向后)。
% 当 delta_mag 接近0.5时，性能最好；接近0时，性能最差。
% 我们设定一个阈值，比如0.25。如果delta_mag大于它，说明次大谱线
% 的能量足够强，可以认为处于“甜点区”。这个值可以根据需求微调。
sweet_spot_threshold = 0.25; 


% 步骤 3: 根据 delta_mag 的值进行智能决策
if  delta_mag > sweet_spot_threshold
    % --- 情况A: 信号落在Rife的甜点区 ---
    % delta_mag 足够大，说明Rife算法的结果是可靠且高精度的。
    % 直接采用这次计算的结果，无需额外操作。
    % disp('决策: 使用 Rife 算法'); % 可选：取消注释以查看决策过程
    f_est = f_rife;
    
else
    % --- 情况B: 信号落在Rife的性能洼地区 ---
    % delta_mag 太小，说明Rife算法的结果可能不准确。
    % 在这种情况下，我们启动为“洼地区”专门优化的M-Rife算法。
    % disp('决策: 使用 M-Rife 算法'); % 可选：取消注释以查看决策过程
    f_est = mrife_algorithm(s, fs);
end

end