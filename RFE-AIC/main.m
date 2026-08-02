function results = main()
% main  QZS-DVA 自适应阻抗控制 — 控制器验证入口
%
%   results = main()
%
%   本阶段唯一目的：快速验证自适应控制器是否有效（单向传递 + simDvaOnly）。
%   耦合宿主、冻结天花板、backbone、自检等留待后续，本轮不接入。

clc; clear; close all;

%% 路径
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
addpath(genpath(fullfile(rootDir, 'lib')));
addpath(genpath(fullfile(thisDir, 'funcs')));
addpath(genpath(fullfile(thisDir, 'sim')));
outDir = fullfile(thisDir, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% 物理与参考模型
% 参考模型参数设置为与物理模型一致，只是刚度去掉了非线性刚度
p     = getQzsDvaParams('isPrintSummary', true);
plant = buildPlant(p);
host  = buildHost(p);        % 仅用于生成宿主传函（力 → 安装点加速度）

%% 时间与采样
sim.fs_ctrl   = 2000;                 % 控制器采样率 Hz
% sim.M 是控制器按较粗时钟更新，RK4仿真对象用更细步长积分，兼顾稳定性与算力
% M 太小：积分误差大、硬非线性可能不稳；太大：更准但更慢，收益通常很快饱和
sim.M         = 8;                    % 每控制周期内对象积分步数
sim.h         = 1/(sim.fs_ctrl*sim.M); % 对象积分步长
sim.T         = 60;                   % 仿真时长 s（验证阶段 60s 足够）
sim.t_discard = 10;                   % 前 10 s 不计入指标

%% 频带选择（二阶带通，相对阶必须为 1）和加速度误差低通滤波器（一阶）
flt.omega_d = 2*pi*2.0;               % 目标主导频率 rad/s
flt.zeta_b  = 0.1;                    % 带通阻尼比，通带越宽、瞬态越快，但带外干扰更多、相位更“黏”
flt.omega_c = 10*flt.omega_d;         % 一阶低通截止频率

%% 控制参数
% 滑模/误差面参数
ctl.w1        = 3;                    % 状态误差面权重：s 中 Be 的权重（Bed + w1*Be + w2*a_f）；原 20 偏大
ctl.lambda    = 0;                  % 反作用力误差权重（与 omega_c 解耦）；过大则零动态下误差难收敛
ctl.w2        = ctl.lambda / flt.omega_c;  % 加速度误差项权重（w2 = lambda/omega_c）
ctl.Ms        = 1 + ctl.lambda;       % = 1+lambda；控制律前的质量归一化因子，Ms越大，同样的u推动 越有力
ctl.K         = 10;                   % 趋近增益：控制律中 -K*s，越大 s→0 越快，但 u 更猛、越易抖振
ctl.Gamma_k   = 0;                  % 等效刚度 Δk 的自适应增益（越大收敛越快，越易漂移）；占位值
ctl.Gamma_c   = 0;                  % 等效阻尼 Δc 的自适应增益
ctl.adapt_on  = true;                 % true=更新 Δk/Δc
ctl.u_sat     = Inf;                  % 控制力饱和限幅；Inf 表示第一版不加饱和
ctl.mode      = 'adaptive';           % 'adaptive'|'frozen'|'passive'(u=0 无控)

%% 激励：力 → 宿主 → 安装点加速度（输出恒为 xdd）
exc.f_d        = 2.0;                 % 主导频率 Hz（宜与 flt.omega_d 一致）
exc.amp        = 0.005;                 % 安装点加速度 RMS 目标值 (m/s^2)——宿主滤波后归一化
exc.rand_ratio = 0;                 % 力端随机功率占比 [0,1]；主线 0.1~0.2
exc.f_lo       = 0.2;                 % 随机带限下界 Hz（设计约束，防 DC 推离 QZS 工作点）
exc.f_hi       = 20.0;                % 随机带限上界 Hz（= 10*f_d，设计约束）
exc.seed       = 1;                   % 随机种子，保证可复现

%% 滤波器离散化
fltSys = designBandpass(flt, sim.fs_ctrl);

%% 激励自检（验证清单 1–3）
tCtrl = (0:1/sim.fs_ctrl:sim.T).';
% 生成的激励是简谐+随机激励下，经过宿主结构后的加速度时程
[xddCheck, excStats] = genExcitation(exc, tCtrl, host, flt);
fprintf('[main] r_inband=%.4f  (应 = rand_ratio=%.3f)\n', ...
    excStats.r_inband, exc.rand_ratio);

figExc = figure('Name', 'Excitation-Check', 'Color', 'w');
tiledlayout(2, 1);
nexttile;
plot(tCtrl, xddCheck); grid on; xlabel('t (s)'); ylabel('xdd (m/s^2)');
title(sprintf('安装点加速度 (RMS=%.3e, r_{inband}=%.3f)', excStats.rms, excStats.r_inband));
nexttile;
plot(excStats.psd_f, excStats.psd_p); grid on;
xlim([0, max(30, exc.f_hi*1.2)]); xlabel('f (Hz)'); ylabel('PSD');
xline(exc.f_d, '--r');
title('激励 PSD（f_d 处应有谱峰）');

%% 公共配置包（传给实验脚本，避免散落全局变量）
cfg.plant       = plant;              % 对象/参考模型物理参数（m,k1,k3,...）
cfg.host        = host;               % 宿主参数（仅用于力→加速度传函）
cfg.sim         = sim;                % 时间与采样（fs_ctrl, M, h, T, t_discard）
cfg.flt         = flt;                % 连续域滤波参数（omega_d, zeta_b, omega_c）
cfg.fltSys      = fltSys;             % 离散滤波器状态空间（bp/lp）
cfg.ctl         = ctl;                % 控制与自适应参数
cfg.exc         = exc;                % 激励默认参数（频率/幅值/种子）
cfg.outDir      = outDir;             % 结果输出目录（.mat/.png）
cfg.tStep       = 30;                 % 自适应场景二：幅值阶跃时刻 s
cfg.ampStepGain = 3;                  % 自适应场景二：阶跃后幅值倍数（×3）
cfg.excStats    = excStats;           % 激励自检统计

results = runAdaptive(cfg);
results.excStats = excStats;

% fprintf('\n[main] 完成。结果目录: %s\n', outDir);
end
