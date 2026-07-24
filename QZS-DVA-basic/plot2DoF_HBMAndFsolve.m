%% plotQzsDvaFreqResponse  一阶 HBM 幅频响应（含正/逆向延拓）
% 基于谐波平衡法代数方程组，用 fsolve 求解并绘制主结构传递率与吸振器位移比。
% 后续可在此脚本末尾追加 ODE 数值积分对比段。

clc
clear
% close all

%% 参数与扫频配置
params = getQzsDvaParams('isPrintSummary', true, 'needDimensionless', true);

R_m = params.R_m;
F_N = params.F_N;
X_st = F_N / params.k0;   % 主结构静位移，用于传递率归一化

% 无量纲频率扫频范围（以 w1Residual 为基准）
Omega_min = 0.05;
Omega_max = 2;
numOmega = 1000;
Omega_vec = linspace(Omega_min, Omega_max, numOmega);

% fsolve 选项：适度容差即可满足幅频曲线绘制
optsFsolve = optimoptions('fsolve', ...
    'Display', 'off', ...
    'FunctionTolerance', 1e-10, ...
    'StepTolerance', 1e-10);

%% 正向扫频（Omega 递增，延拓捕捉上分支）
Z_init = zeros(4, 1);
[Omega_fwd, Z_fwd, isOk_fwd] = sweepHbmBranch(Omega_vec, Z_init, optsFsolve, params);

%% 逆向扫频（Omega 递减，延拓捕捉下分支与跳跃区）
Z_init = Z_fwd(end, :).';   % 从正向末端解出发，便于落入不同吸引子
[Omega_bwd, Z_bwd, isOk_bwd] = sweepHbmBranch(flipud(Omega_vec), Z_init, optsFsolve, params);

%% 物理量还原  trans_fwd: 主结构位移传递率 ratio_fwd: 吸振器绝对位移 / 主结构位移
[f_fwd, trans_fwd, localAcc2Disp_fwd, localDispRatio_fwd, B_eq_fwd] = ...
    restorePhysicalAmplitudes(Omega_fwd, Z_fwd, params, X_st);
[f_bwd, trans_bwd, localAcc2Disp_bwd, localDispRatio_bwd, B_eq_bwd] = ...
    restorePhysicalAmplitudes(Omega_bwd, Z_bwd, params, X_st);

%% 绘图
% 图 1：主结构位移传递率 A_x0 / X_st
fig1 = figure('Name', '主结构位移传递率', 'Color', 'w');
hold on;
plot(f_fwd, 20*log10(trans_fwd), 'b-', 'LineWidth', 1.2, 'DisplayName', '正向扫频');
plot(f_bwd, 20*log10(trans_bwd), 'r--', 'LineWidth', 1.2, 'DisplayName', '逆向扫频');
grid on;
xlabel('频率 f (Hz)');
ylabel('主结构位移传递率 A_{x0} / X_{st}');
title('主结构位移传递率');
legend('Location', 'best');

% 图 2：宿主结构加速度到局部位移传递率 A_y / A_ddot_x0bar
fig2 = figure('Name', '宿主加速度到局部位移传递率', 'Color', 'w');
hold on;
plot(f_fwd, 20*log10(localAcc2Disp_fwd), 'b-', 'LineWidth', 1.2, 'DisplayName', '正向扫频');
plot(f_bwd, 20*log10(localAcc2Disp_bwd), 'r--', 'LineWidth', 1.2, 'DisplayName', '逆向扫频');
grid on;
xlabel('频率 f (Hz)');
ylabel('20log_{10}(|A_y|/|A_{a0}|) (dB)');
title('宿主结构加速度到局部位移传递率');
legend('Location', 'best');

%% ===== 局部函数 =====

function [Omega_out, Z_out, isOk_out] = sweepHbmBranch(Omega_vec, Z_init, optsFsolve, params)
% 沿给定 Omega 序列做延拓扫频，上一步解作为下一步初值

numOmega = numel(Omega_vec);
Z_out = nan(numOmega, 4);
isOk_out = false(numOmega, 1);
Z_curr = Z_init;

for iOmega = 1:numOmega
    Omega = Omega_vec(iOmega);
    hbmFun = @(Z) hbmResidual(Z, Omega, params);
    [Z_sol, Fval, exitflag] = fsolve(hbmFun, Z_curr, optsFsolve);

    exitflag_out(iOmega) = exitflag;
    resnorm_out(iOmega) = norm(Fval);

    if exitflag > 0 && norm(Fval) < 1e-7
        Z_out(iOmega, :) = Z_sol.';
        isOk_out(iOmega) = true;
        Z_curr = Z_sol;
    end
    % 若求解失败则保留 NaN，但不重置 Z_curr，尽量维持延拓链
end

Omega_out = Omega_vec;
end

function F = hbmResidual(Z, Omega, params)
% 一阶 HBM 四方程残差，Z = [u1; v1; u2; v2]

u1 = Z(1);
v1 = Z(2);
u2 = Z(3);
v2 = Z(4);

mu = params.mu;
lambda = params.lambda;
alpha_n3 = params.alpha_n3;
alpha_0 = params.alpha0;
xi0 = params.xi0;
xi1 = params.xi1;
X_force = params.X_force;

A1_sq = u1^2 + v1^2;
nl_u1 = (3 / 4) * alpha_n3 * u1 * A1_sq;
nl_v1 = (3 / 4) * alpha_n3 * v1 * A1_sq;

F = zeros(4, 1);
F(1) = (1 - Omega^2) * u1 - Omega^2 * u2 + 2 * xi1 * Omega * v1 + nl_u1;
F(2) = (1 - Omega^2) * v1 - Omega^2 * v2 - 2 * xi1 * Omega * u1 + nl_v1;
F(3) = ( mu * alpha_0 - Omega^2) * u2 - mu * u1 ...
    + 2 * xi0 * lambda * Omega * v2 - 2 * mu * xi1 * Omega * v1 ...
    - (3 / 4) * mu * alpha_n3 * u1 * A1_sq;
F(4) = ( mu * alpha_0 - Omega^2) * v2 - mu * v1 ...
    - 2 * xi0 * lambda * Omega * u2 + 2 * mu * xi1 * Omega * u1 ...
    - (3 / 4) * mu * alpha_n3 * v1 * A1_sq - X_force;
end

function [f_Hz, transRatio, localAcc2Disp, localDispRatio, B_eq] = ...
    restorePhysicalAmplitudes(Omega_vec, Z_mat, params, X_st)
% transRatio      = A_x0 / X_st，主结构位移传递率
% localAcc2Disp   = A_y / A_ddot_x0bar，无量纲宿主加速度到局部位移传递率
% localDispRatio  = A_y / A_x0，局部位移 / 主结构位移
% B_eq            = A_ddot_x0bar = Omega^2 * A_x0bar，等效宿主加速度输入幅值

R_m = params.R_m;
w1Residual = params.w1Residual;
numPts = size(Z_mat, 1);

f_Hz = nan(numPts, 1);
transRatio = nan(numPts, 1);
localAcc2Disp = nan(numPts, 1);
localDispRatio = nan(numPts, 1);
B_eq = nan(numPts, 1);

for iPt = 1:numPts
    u1 = Z_mat(iPt, 1);
    v1 = Z_mat(iPt, 2);
    u2 = Z_mat(iPt, 3);
    v2 = Z_mat(iPt, 4);

    if any(isnan([u1, v1, u2, v2]))
        continue;
    end

    Omega = Omega_vec(iPt);
    f_Hz(iPt) = Omega * w1Residual / (2 * pi);

    % 无量纲幅值
    A_y_bar  = hypot(u1, v1);     % 局部相对位移 y 的无量纲幅值
    A_x0_bar = hypot(u2, v2);     % 主结构位移 x0/R 的无量纲幅值

    % 物理主结构位移幅值
    A_x0 = R_m * A_x0_bar;

    % 主结构位移传递率
    transRatio(iPt) = A_x0 / X_st;

    % 等效宿主结构无量纲加速度输入幅值
    % 若 x0bar = u2 cos(Omega tau) + v2 sin(Omega tau)
    % 则 |ddot{x0bar}| = Omega^2 * sqrt(u2^2 + v2^2)
    B_eq(iPt) = Omega^2 * A_x0_bar;

    % 局部位移 / 主结构位移
    if A_x0_bar > 1e-12
        localDispRatio(iPt) = A_y_bar / A_x0_bar;
    end

    % 宿主结构加速度 -> 局部位移传递率  
    % 为什么需要除以 w1Residual^2？ 因为在无量纲公式推导下的时间尺度单位是τ，不是t，物理量纲下的传递率是t，所以需要除以w1Residual^2转换到物理量纲。
    if B_eq(iPt) > 1e-12
        localAcc2Disp(iPt) = A_y_bar / B_eq(iPt) / w1Residual^2;   % 单位: m / (m/s^2) = s^2/m
    end
end
end