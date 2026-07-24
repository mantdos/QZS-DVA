%% plotSdoFQzsDvaLocalFreqResponse  单自由度 QZS-DVA 局部方程一阶 HBM 幅频响应
% 基于局部无量纲方程
%   y'' + 2*xi1*y' + y + alpha_n3*y^3 = -x0_bar''
% 在宿主结构加速度简谐输入下，用一阶谐波平衡法导出幅值方程，并化为
% 关于 s = A^2 的三次代数方程，用 roots 求全部正实根以检查参数合理性。
%
% 说明:
%   B — 宿主结构加速度输入的无量纲幅值（不是外部激励力幅值）
%   A — 吸振器局部相对位移 y 的无量纲幅值
%   H = A/B — 宿主加速度到局部位移的幅值传递率（无量纲 FRF 模）
%   同一频率点若出现多个正实根，表示一阶 HBM 下可能存在多值幅频响应
%   （跳跃、滞回等非线性现象的典型征兆）

clc
clear
% close all

%% 参数与扫频配置
params = getQzsDvaParams('needDimensionless', true);

xi1       = params.xi1;
alpha_n3  = params.alpha_n3;
w1Residual = params.w1Residual;

% 无量纲频率扫频范围（与 plotQzsDvaFreqResponse.m 保持一致）
Omega_min = 0.05;
Omega_max = 2;
numOmega  = 1000;
Omega_vec = linspace(Omega_min, Omega_max, numOmega);

% 宿主结构无量纲加速度输入幅值向量（可按需修改）
A_a0_vec = [1e-2, 5e-2, 1e-1, 5e-1];

% 三次方程求根时的虚部容差
imagTol = 1e-8;

beta = (3 / 4) * alpha_n3;

%% 线性系统对照（alpha_n3 = 0 时的解析传递率，与 B 无关）
denom_linear = sqrt((1 - Omega_vec.^2).^2 + (2 * xi1 .* Omega_vec).^2);
H_linear = 1 ./ denom_linear / w1Residual^2;   % H = A/B = 1 / |1 - Omega^2 + j*2*xi1*Omega|

%% 非线性 HBM：逐频率、逐输入幅值求三次方程全部正实根
% 存储结构: 每个 B 对应 Omega/A 数组（含多值分支）
hbmResults = struct('B', {}, 'Omega', {}, 'A', {}, 'numRoots', {});

fprintf('\n===== 单自由度局部 HBM 幅频扫描 =====\n');
fprintf('xi1 = %.4f, alpha_n3 = %.4e, w1Residual = %.4f rad/s\n', ...
    xi1, alpha_n3, w1Residual);

for iB = 1:numel(A_a0_vec)
    A_a0 = A_a0_vec(iB);

    Omega_all = [];
    A_y_all     = [];
    numRoots_perOmega = zeros(numOmega, 1);

    for iOm = 1:numOmega
        Omega = Omega_vec(iOm);

        % 幅值方程 B^2 = s * [(1-Omega^2 + beta*s)^2 + (2*xi1*Omega)^2]
        % 整理为 a3*s^3 + a2*s^2 + a1*s + a0 = 0
        a3 = beta^2;
        a2 = 2 * (1 - Omega^2) * beta;
        a1 = (1 - Omega^2)^2 + (2 * xi1 * Omega)^2;
        a0 = -A_a0^2;

        s_roots = roots([a3, a2, a1, a0]);

        % 只保留正实根: s = A^2 > 0，非正实根在数学上没有含义
        isPosReal = abs(imag(s_roots)) < imagTol & real(s_roots) > 0;
        s_valid = real(s_roots(isPosReal));
        nValid = numel(s_valid);
        numRoots_perOmega(iOm) = nValid;

        if nValid > 0
            A_y_valid = sqrt(s_valid); % 每个Omega下的正实根A值
            Omega_all = [Omega_all; repmat(Omega, nValid, 1)]; % 如果存在多个正实根，则填入多个Omega对应多个A值
            A_y_all     = [A_y_all;     A_y_valid];                  
        end
    end

    hbmResults(iB).A_a0         = A_a0;
    hbmResults(iB).Omega     = Omega_all;
    hbmResults(iB).A_y         = A_y_all;
    hbmResults(iB).numRoots  = numRoots_perOmega;

    % 命令行输出：出现多个正实根的频率范围
    multiMask = numRoots_perOmega > 1;
    fprintf('\n--- A_a0 = %.4g ---\n', A_a0);
    if ~any(multiMask)
        fprintf('  全频段仅存在单值 HBM 解（每个 Omega 至多 1 个正实根）。\n');
    else
        ranges = findContiguousRanges(Omega_vec, multiMask);
        fprintf('  存在多个正实根的频率区间 (Omega):\n');
        for iR = 1:size(ranges, 1)
            omLo = ranges(iR, 1);
            omHi = ranges(iR, 2);
            fLo  = omLo * w1Residual / (2 * pi);
            fHi  = omHi * w1Residual / (2 * pi);
            fprintf('    Omega in [%.4f, %.4f]  (f in [%.2f, %.2f] Hz)\n', ...
                omLo, omHi, fLo, fHi);
        end
    end
end

%% 图 1：无量纲传递率曲线  Omega — H = A/B
fig1 = figure('Name', '局部 HBM 无量纲传递率', 'Color', 'w');
hold on;
colors = lines(numel(A_a0_vec));

for iB = 1:numel(A_a0_vec)
    res = hbmResults(iB);
    H_hbm = res.A_y / res.A_a0 / w1Residual^2;   % 宿主加速度 → 局部位移传递率

    % 非线性 HBM 全部正实根（散点，保留多值分支）
    if ~isempty(res.Omega)
        scatter(res.Omega, 20 * log10(H_hbm), 8, colors(iB, :), 'filled', ...
            'DisplayName', sprintf('HBM, A_a0=%.4g', res.A_a0));
    end
end

% 线性对照（与 B 无关，只画一条）
plot(Omega_vec, 20 * log10(H_linear), 'k--', 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

grid on;
xlabel('\Omega');
ylabel('位移/加速度传递率 |A/B| (dB)');
title('单自由度 QZS-DVA 局部方程 — 一阶 HBM 无量纲传递率');
legend('Location', 'best');

%% 图 2：物理频率传递率曲线  f (Hz) — |y|/|a_0|
% 物理幅值: y = A*R_m, a_0 = B*w1Residual^2*R_m  =>  |y|/|a_0| = A/(B*w1Residual^2)
fig2 = figure('Name', '局部 HBM 物理传递率', 'Color', 'w');
hold on;

for iB = 1:numel(A_a0_vec)
    res = hbmResults(iB);
    % 为什么需要除以 w1Residual^2？ 因为在无量纲公式推导下的时间尺度单位是τ，不是t，物理量纲下的传递率是t，所以需要除以w1Residual^2转换到物理量纲。
    H_phys = res.A_y ./ (res.A_a0 * w1Residual^2);   % 单位: m / (m/s^2) = s^2/m

    f_Hz = res.Omega * w1Residual / (2 * pi);

    if ~isempty(res.Omega)
        scatter(f_Hz, 20 * log10(H_phys), 8, colors(iB, :), 'filled', ...
            'DisplayName', sprintf('HBM, A_a0=%.4g', res.A_a0));
    end
end

f_linear = Omega_vec * w1Residual / (2 * pi);
plot(f_linear, 20 * log10(H_linear), 'k--', 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

grid on;
xlabel('频率 f (Hz)');
ylabel('位移/加速度传递率 |y|/|a_0| (dB)');
title('单自由度 QZS-DVA 局部方程 — 一阶 HBM 物理传递率');
legend('Location', 'best');

%% ===== 局部函数 =====

function ranges = findContiguousRanges(xVec, mask)
% 在逻辑向量 mask 为 true 的区段上，返回 xVec 的连续区间 [xLo, xHi]
% 用于报告多值响应出现的频率范围

if ~any(mask)
    ranges = zeros(0, 2);
    return;
end

mask = mask(:);
edges = diff([false; mask; false]);
iStart = find(edges == 1);
iEnd   = find(edges == -1) - 1;

nSeg = numel(iStart);
ranges = zeros(nSeg, 2);
for k = 1:nSeg
    ranges(k, 1) = xVec(iStart(k));
    ranges(k, 2) = xVec(iEnd(k));
end
end
