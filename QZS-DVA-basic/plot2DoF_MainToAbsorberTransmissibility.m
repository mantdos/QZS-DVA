%% plot2DoF_MainToAbsorberTransmissibility  主系统→吸振器传递率（线性二自由度）
% 宿主结构到吸振器质量的频率响应:
%   H(s) = (c*s + k) / ((m1 + ma)*s^2 + c*s + k)
% 其中 ma 为可调附加质量，s = j*omega。
%
% 物理参数取自图 1-4 二自由度系统算例。

clc
clear
% close all

%% 物理参数（图 1-4）
m1 = 61e-3;          % 宿主结构等效质量 (kg)
k  = 2.3727e4;       % 刚度 (N/m)
c  = 2.5;            % 阻尼系数 (Ns/m)

%% 扫频与 ma 取值
f_min   = 50;        % Hz
f_max   = 150;       % Hz
numFreq = 500;
f_Hz    = linspace(f_min, f_max, numFreq);
omega   = 2 * pi * f_Hz;

% 附加质量 ma (kg)，可按需修改
ma_vec = [-0.02,-0.01,0, 0.01,  0.03,  0.05];

%% 逐 ma 计算传递率幅值 |H(j*omega)|
fig = figure('Name', '主系统→吸振器传递率', 'Color', 'w');
hold on;

colors = lines(numel(ma_vec));

for iMa = 1:numel(ma_vec)
    ma = ma_vec(iMa);
    m_eq = m1 + ma;

    % H(j*omega) = (c*j*omega + k) / (-m_eq*omega^2 + c*j*omega + k)
    num = c * 1i * omega + k;
    den = -m_eq * omega.^2 + c * 1i * omega + k;
    H_mag = abs(num ./ den);

    plot(f_Hz, 20 * log10(H_mag), 'LineWidth', 1.5, ...
        'Color', colors(iMa, :), ...
        'DisplayName', sprintf('m_a = %.3f kg', ma));
end

grid on;
xlabel('频率 f (Hz)');
ylabel('传递率 |H(j\omega)| (dB)');
title('主系统 \rightarrow 吸振器 频率响应（不同 m_a）');
xlim([f_min, f_max]);
legend('Location', 'best');

fprintf('\n===== 主系统→吸振器线性传递率 =====\n');
fprintf('m_1 = %.4f kg, k = %.4f N/m, c = %.4f Ns/m\n', m1, k, c);
fprintf('频率范围: %.1f – %.1f Hz\n', f_min, f_max);
for iMa = 1:numel(ma_vec)
    ma = ma_vec(iMa);
    fn = sqrt(k / (m1 + ma)) / (2 * pi);
    fprintf('  m_a = %.3f kg  ->  固有频率 f_n = %.2f Hz\n', ma, fn);
end
