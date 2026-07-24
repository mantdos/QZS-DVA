%% plot2DoF_HBMAndRoots  二自由度 QZS-DVA 一阶 HBM 幅频响应（解析求根法）
% 基于一阶谐波平衡法（HBM）代数消元，将四方程组化为关于幅值平方 s 的一元
% 三次方程，用 MATLAB roots 一次性求出每个频率下的全部正实根，从而完整
% 展现非线性多值分支（跳跃、滞回等现象），无需 fsolve 数值延拓。
%
% 三次方程（对每个 Omega）:
%   P2*s^3 + P1*s^2 + P0*s - Omega^4 * X_force^2 = 0
% 其中 s = A1^2，A1 为吸振器局部相对位移 y 的无量纲幅值。
%
% 对每个合法正实根 s:
%   K   = k0 + k1*s          等效线性刚度
%   A2  = sqrt(s)*sqrt(K^2+C^2) / Omega^2   宿主结构无量纲位移幅值
%   transRatio = A2 / X_force               物理位移传递率
%
% 参数、扫频范围及绘图格式与 plot2DoF_HBMAndFsolve.m 保持一致。

clc
clear
% close all

%% 参数与扫频配置
params = getQzsDvaParams('isPrintSummary', true, 'needDimensionless', true);

R_m        = params.R_m;
F_N        = params.F_N;
X_st       = F_N / params.k0;          % 主结构静位移，与 fsolve 脚本一致
w1Residual = params.w1Residual;

mu       = params.mu;
lambda   = params.lambda;
alpha_n3 = params.alpha_n3;
alpha_0  = params.alpha0;
xi0      = params.xi0;
xi1      = params.xi1;
X_force  = params.X_force;


% 无量纲频率扫频范围（以 w1Residual 为基准，与 plot2DoF_HBMAndFsolve.m 一致）
Omega_min = 0.05;
Omega_max = 2;
numOmega  = 200;
Omega_vec = linspace(Omega_min, Omega_max, numOmega);

% 求根时判定“实根”的虚部容差
imagTol = 1e-8;

% 是否运行 ODE45 数值仿真验证（初值延续法）
runNumericalSim = false;

%% 逐频率解析求根，收集全部正实根分支
Omega_all       = [];
A_0_all         = [];   % 宿主结构无量纲位移幅值
transRatio_all  = [];   % 激励力到宿主结构位移传递率 A_0 / X_force
acc2LocalDisp_all = []; % 宿主结构加速度到局部位移传递率 A_y / A_ddot_x0bar
numRoots_perOmega = zeros(numOmega, 1);

fprintf('\n===== 二自由度 QZS-DVA HBM 解析求根扫频 =====\n');
fprintf('mu=%.4f, lambda=%.4f, alpha_n3=%.4e, X_force=%.4f\n', ...
    mu, lambda, alpha_n3, X_force);

for iOm = 1:numOmega
    Omega = Omega_vec(iOm);

    % --- 基础变量 ---
    k0 = 1 - Omega^2;
    k1 = (3 / 4) * alpha_n3;
    C  = 2 * xi1 * Omega;

    H22 = mu * alpha_0 - (1 + mu) * Omega^2;
    H21 = -2 * xi0 * lambda * Omega;

    H22_sq = H22^2 + H21^2;

    % --- 三次多项式系数 ---
    P2 = k1^2 * H22_sq;
    P1 = 2 * k0 * k1 * H22_sq - 2 * mu * Omega^4 * k1 * H22;
    P0 = (k0^2 + C^2) * H22_sq ...
        - 2 * mu * Omega^4 * k0 * H22 ...
        - 2 * mu * Omega^4 * C * H21 ...
        + mu^2 * Omega^8;

    % 常数项: -Omega^4 * X_force^2
    polyCoeffs = [P2, P1, P0, -Omega^4 * X_force^2];
    s_roots = roots(polyCoeffs);

    % --- 筛选正实根 s = A1^2 > 0 ---
    isPosReal = abs(imag(s_roots)) < imagTol & real(s_roots) > 0;
    s_valid = real(s_roots(isPosReal));
    nValid = numel(s_valid);
    numRoots_perOmega(iOm) = nValid;

    if nValid == 0
        continue;
    end

    % 对每个合法根计算宿主结构幅值及传递率
    for iRoot = 1:nValid
        s = s_valid(iRoot);

        K  = k0 + k1 * s;
        A_y = sqrt(s); % sqrt{u_1^2 + v_1^2}，无量纲局部位移幅值

        % 宿主结构无量纲位移幅值
        if Omega < 1e-12
            continue;
        end
        A_0 = A_y * sqrt(K^2 + C^2) / Omega^2; % 无量纲主结构位移幅值

        if A_0 <= 0
            continue;
        end

        % 物理位移传递率（按 HBM 消元结果定义）
        transRatio = A_0 / X_force;

        % 与 fsolve 脚本一致的附加物理量（用于第二幅图）
        B_eq    = Omega^2 * A_0; % 等效无量纲宿主加速度输入幅值
        if B_eq > 1e-12
            % 为什么需要除以 w1Residual^2？ 因为在无量纲公式推导下的时间尺度单位是τ，不是t，物理量纲下的传递率是t，所以需要除以w1Residual^2转换到物理量纲。
            acc2LocalDisp = A_y / B_eq / w1Residual^2; % 单位: m / (m/s^2) = s^2/m
        else
            acc2LocalDisp = NaN;
        end

        Omega_all         = [Omega_all;         Omega]; %#ok<AGROW>
        A_0_all            = [A_0_all;            A_0];    %#ok<AGROW>
        transRatio_all    = [transRatio_all;    transRatio]; %#ok<AGROW>
        acc2LocalDisp_all = [acc2LocalDisp_all; acc2LocalDisp]; %#ok<AGROW>
    end
end

%% 命令行报告：出现多个正实根的频率区间
multiMask = numRoots_perOmega > 1;
if ~any(multiMask)
    fprintf('全频段每个 Omega 至多 1 个正实根（HBM 单值解）。\n');
else
    ranges = findContiguousRanges(Omega_vec, multiMask);
    fprintf('存在多个正实根（多值/跳跃区）的频率区间:\n');
    for iR = 1:size(ranges, 1)
        omLo = ranges(iR, 1);
        omHi = ranges(iR, 2);
        fLo  = omLo * w1Residual / (2 * pi);
        fHi  = omHi * w1Residual / (2 * pi);
        fprintf('  Omega in [%.4f, %.4f]  (f in [%.2f, %.2f] Hz)\n', ...
            omLo, omHi, fLo, fHi);
    end
end

%% ODE45 数值仿真（初值延续法，验证解析求根）
if runNumericalSim
    odeFuncs = FuncsQzsOde();

    numPeriods     = 100; % 每个频率点仿真的激励周期数（保证瞬态衰减）
    extractPeriods = 10; % 提取稳态周期数（保证稳态）
    odeOpts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

    A_y_fwd  = zeros(numOmega, 1);
    A_0_fwd  = zeros(numOmega, 1);

    fprintf('\n===== ODE45 数值仿真扫频（正向） =====\n');

    % 选定频率点（Hz）→ 无量纲 Omega，并在扫频网格上取最近索引
    f_Hz_plot = [2.17, 2.3];
    Omega_plot_vec = f_Hz_plot * (2 * pi) / w1Residual;
    Omega_plot_idx = zeros(size(f_Hz_plot));
    for kPlot = 1:numel(f_Hz_plot)
        [~, Omega_plot_idx(kPlot)] = min(abs(Omega_vec - Omega_plot_vec(kPlot)));
    end

    % --- 正向扫频：Omega 递增，零初值出发 ---
    Z_init = [0; 0; 0; 0];
    for iOm = 1:numOmega
        Omega = Omega_vec(iOm);
        T_period = 2 * pi / Omega; % 当前频率简谐激励完成一个周期的时间
        tau_span = [0, numPeriods * T_period]; % 每个频率下的总仿真时间
        
        % t_sol：求解器自动选取的时间点，Z_sol：对应的状态
        [t_sol, Z_sol] = ode45(@(t, Z) odeFuncs.ode3rd(t, Z, Omega, params), ...
            tau_span, Z_init, odeOpts);

        tau_extract = (numPeriods - extractPeriods) * T_period; % 提取稳态周期开始的时间点
        idxSteady = t_sol >= tau_extract; % 提取稳态周期的时间点索引
        y_steady  = Z_sol(idxSteady, 1);
        x0_steady = Z_sol(idxSteady, 3);
        
        % 通过最大值和最小值计算稳态幅值
        A_y_fwd(iOm) = (max(y_steady)  - min(y_steady))  / 2;
        A_0_fwd(iOm) = (max(x0_steady) - min(x0_steady)) / 2;

        % 在选定的两个频率点绘制稳态时域响应
        for kPlot = 1:numel(Omega_plot_idx)
            if iOm ~= Omega_plot_idx(kPlot)
                continue;
            end
            Om_plot = Omega_vec(iOm);
            f_plot  = Om_plot * w1Residual / (2 * pi);
            fprintf('展示扫频中间结果：Omega = %.4f, frequency = %.2f Hz\n', Om_plot, f_plot);
            figure(100 + kPlot);
            clf;
            ax1 = subplot(2, 1, 1);
            plot(ax1, t_sol(idxSteady), y_steady);
            ylabel(ax1, '局部位移 y');
            grid(ax1, 'on');
            ax2 = subplot(2, 1, 2);
            plot(ax2, t_sol(idxSteady), x0_steady);
            ylabel(ax2, '主结构位移 x0');
            xlabel(ax2, '时间 \tau');
            grid(ax2, 'on');
            title(sprintf('扫频中间结果：Omega = %.4f，frequency = %.2f Hz', Om_plot, f_plot));
        end

        Z_init = Z_sol(end, :).';
    end

    fprintf('===== ODE45 数值仿真扫频（逆向） =====\n');

    % --- 逆向扫频：Omega 递减，承接正向末端状态 ---
    Omega_bwd_vec = flipud(Omega_vec(:));
    A_y_bwd_raw = zeros(numOmega, 1);
    A_0_bwd_raw = zeros(numOmega, 1);

    for iOm = 1:numOmega
        Omega = Omega_bwd_vec(iOm);
        T_period = 2 * pi / Omega;
        tau_span = [0, numPeriods * T_period];

        [t_sol, Z_sol] = ode45(@(t, Z) odeFuncs.ode3rd(t, Z, Omega, params), ...
            tau_span, Z_init, odeOpts);

        tau_extract = (numPeriods - extractPeriods) * T_period;
        idxSteady = t_sol >= tau_extract;
        y_steady  = Z_sol(idxSteady, 1);
        x0_steady = Z_sol(idxSteady, 3);

        A_y_bwd_raw(iOm) = (max(y_steady)  - min(y_steady))  / 2;
        A_0_bwd_raw(iOm) = (max(x0_steady) - min(x0_steady)) / 2;

        Z_init = Z_sol(end, :).';
    end

    % 逆向结果反转，与 Omega_vec 坐标对齐
    A_y_bwd = flipud(A_y_bwd_raw);
    A_0_bwd = flipud(A_0_bwd_raw);

    % --- 传递率指标 ---
    % 激励力到宿主结构位移传递率
    transRatio_fwd     = A_0_fwd / X_force;
    transRatio_bwd     = A_0_bwd / X_force;
    % 等效无量纲宿主加速度输入幅值
    B_eq_fwd           = (Omega_vec').^2 .* A_0_fwd;
    B_eq_bwd           = (Omega_vec').^2 .* A_0_bwd;
    % 宿主结构加速度到局部位移传递率
    acc2LocalDisp_fwd  = A_y_fwd ./ B_eq_fwd / w1Residual^2;
    acc2LocalDisp_bwd  = A_y_bwd ./ B_eq_bwd / w1Residual^2;

    f_Hz_ode = Omega_vec * w1Residual / (2 * pi);
end

%% 物理频率（Hz）
f_Hz_all = Omega_all * w1Residual / (2 * pi);

%% 绘图
% 图 1：主结构位移传递率 A2 / X_force — 散点展示全部正实根分支
fig1 = figure('Name', '主结构位移传递率 (HBM 解析求根)', 'Color', 'w');
hold on;
validMask1 = transRatio_all > 0 & isfinite(transRatio_all);
scatter(f_Hz_all(validMask1), 20 * log10(transRatio_all(validMask1)), ...
    8, [0.2, 0.4, 0.9], 'filled', 'DisplayName', 'HBM 全部正实根');
grid on;
xlabel('频率 f (Hz)');
ylabel('主结构位移传递率 A_2 / X_{force} (dB)');
title('主结构位移传递率 — 一阶 HBM 解析求根（含多值分支）');
if runNumericalSim
    plot(f_Hz_ode, 20 * log10(transRatio_fwd), 'k-', 'LineWidth', 1.5, ...
        'DisplayName', '数值仿真-正向');
    plot(f_Hz_ode, 20 * log10(transRatio_bwd), 'g--', 'LineWidth', 1.5, ...
        'DisplayName', '数值仿真-逆向');
end
legend('Location', 'best');

% 图 2：宿主结构加速度到局部位移传递率
fig2 = figure('Name', '宿主加速度到局部位移传递率 (HBM 解析求根)', 'Color', 'w');
hold on;
validMask2 = acc2LocalDisp_all > 0 & isfinite(acc2LocalDisp_all);
scatter(f_Hz_all(validMask2), 20 * log10(acc2LocalDisp_all(validMask2)), ...
    8, [0.85, 0.3, 0.2], 'filled', 'DisplayName', 'HBM 全部正实根');
grid on;
xlabel('频率 f (Hz)');
ylabel('20log_{10}(|A_y|/|A_{a0}|) (dB)');
title('宿主结构加速度到局部位移传递率 — 一阶 HBM 解析求根');
if runNumericalSim
    plot(f_Hz_ode, 20 * log10(acc2LocalDisp_fwd), 'k-', 'LineWidth', 1.5, ...
        'DisplayName', '数值仿真-正向');
    plot(f_Hz_ode, 20 * log10(acc2LocalDisp_bwd), 'g--', 'LineWidth', 1.5, ...
        'DisplayName', '数值仿真-逆向');
end
legend('Location', 'best');

analyze2DoF_DynamicBehavior(2.18*2*pi/w1Residual, 1, params);














%% ===== 局部函数 =====

function ranges = findContiguousRanges(xVec, mask)
% 在逻辑向量 mask 为 true 的区段上，返回 xVec 的连续区间 [xLo, xHi]

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
