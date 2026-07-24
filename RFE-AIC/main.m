%% main.m  QZS-DVA 自适应阻抗控制仿真（第一版入口）
%
% 执行顺序（与 readme 第 11 节一致）:
%   1) 激励 + 包络诊断
%   2) 无控幅值-反共振漂移
%   3) runFrozenSweep 天花板测试
%   4) runAdaptive 三场景（天花板通过后再跑）
%
% 用法:
%   main                  % 默认跑 step1~2，并询问是否继续
%   main('steps', 1:3)    % 指定步骤
%   main('steps', 1:4, 'quick', true)  % 快速模式（少种子/短时/不扫频）

function main(options)
arguments
    options.steps (1,:) double = [1]
    options.quick (1,1) logical = false
end

%% 路径
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
addpath(genpath(fullfile(rootDir, 'lib')));
addpath(genpath(fullfile(thisDir, 'funcs')));
addpath(genpath(fullfile(thisDir, 'sim')));
outDir = fullfile(thisDir, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% 物理与参考模型
p = getQzsDvaParams('isPrintSummary', true);
plant = buildPlant(p);

%% 时间与采样
sim.fs_ctrl   = 2000;                 % 控制器采样率 Hz
% sim.M 是控制器按较粗时钟更新，RK4仿真对象用更细步长积分，兼顾稳定性与算力
% M 太小：积分误差大、硬非线性可能不稳；太大：更准但更慢，收益通常很快饱和
sim.M         = 8;                    % 每控制周期内对象积分步数
sim.h         = 1/(sim.fs_ctrl*sim.M); % 对象积分步长
sim.T         = 120;                  % 仿真时长 s
sim.t_discard = 10;                   % 前 10 s 不计入指标

%% 频带选择（二阶带通，相对阶必须为 1）和加速度误差低通滤波器（一阶）
flt.omega_d = 2*pi*2.0;               % 目标主导频率 rad/s
flt.zeta_b  = 0.1;                    % 带通阻尼比，通带越宽、瞬态越快，但带外干扰更多、相位更“钝”
flt.omega_c = 10*flt.omega_d;         % 一阶低通截止频率

%% 控制参数
% 滑模/误差面参数
ctl.c1       = 20;                    % 状态误差面系数：s 中 Be 的权重（Bed + c1*Be + c2*a_f）
ctl.lambda   = 1.0;                   % 反作用力误差权重（与 omega_c 解耦）；过大则零动态下误差难收敛
ctl.c2       = ctl.lambda / flt.omega_c;  % 加速度误差项系数（c2 = lambda/omega_c）
ctl.Ms       = 1 + ctl.c2*flt.omega_c;    % = 1+lambda；控制律前的质量归一化因子，Ms越大，同样的u推动 越有力
ctl.K        = 50;                    % 趋近增益：控制律中 -K*s，越大 s→0 越快，但 u 更猛、越易抖振
ctl.Gamma_k  = 1e3;                   % 等效刚度 Δk 的自适应增益（越大收敛越快，越易漂移）
ctl.Gamma_c  = 1e2;                   % 等效阻尼 Δc 的自适应增益
ctl.adapt_on = true;                  % true=更新 Δk/Δc；false=冻结，使用 dk_frozen/dc_frozen
ctl.dk_frozen = 0;                    % 冻结模式下的等效刚度增量（天花板测试时填理论值）
ctl.dc_frozen = 0;                    % 冻结模式下的等效阻尼增量
ctl.u_sat    = Inf;                   % 控制力饱和限幅；Inf 表示第一版不加饱和
ctl.mode     = 'adaptive';            % 'adaptive'|'frozen'|'passive'(u=0 无控)

%% 激励（基座加速度 xdd；amp 按总 RMS 标定）
exc.f_d        = 2.0;                 % 主导频率 Hz（宜与 flt.omega_d 一致）
exc.bw         = 0.4;                 % 随机分量窄带带宽 Hz；与 rand_ratio 同为 0 时退化为纯简谐
exc.rand_ratio = 0.3;                 % 随机分量功率占比 [0,1]；0=纯简谐，越大包络起伏越明显
exc.amp        = 1.0;                 % 总 RMS 幅值标度（m/s^2）
exc.seed       = 1;                   % 随机种子，保证可复现

%% 快速模式缩短时长 / 减少工况
if options.quick
    sim.T = 40;
    sim.t_discard = 8;
    nSeed = 2;
    ampList = logspace(-1.5, -0.5, 5);
    doSweepId = false;
    fprintf('[main] 快速模式: T=%.0fs, nSeed=%d, 不扫频识别\n', sim.T, nSeed);
else
    nSeed = 1;
    ampList = logspace(-2, 0, 10);
    doSweepId = true;
end

%% 滤波器离散化
fltSys = designBandpass(flt, sim.fs_ctrl);
tau_bp = 1 / (flt.zeta_b * flt.omega_d);   % ≈ 2/(2*zeta*wd)
fprintf('[main] 带通时间常数 ≈ %.3f s\n', tau_bp);
fprintf('[main] 线性固有频率 f1 = %.4f Hz\n', plant.f1_Hz);

%% 公共配置包（传给各实验脚本，避免散落全局变量）
cfg.plant  = plant;                   % 对象/参考模型物理参数（m,k1,k3,...）
cfg.sim    = sim;                     % 时间与采样（fs_ctrl, M, h, T, t_discard）
cfg.flt    = flt;                     % 连续域滤波参数（omega_d, zeta_b, omega_c）
cfg.fltSys = fltSys;                  % 离散滤波器状态空间（bp/lp）
cfg.ctl    = ctl;                     % 控制与自适应参数
cfg.exc    = exc;                     % 激励默认参数（频率/带宽/幅值/种子）
cfg.outDir = outDir;                  % 结果输出目录（.mat/.fig/.png）
cfg.ampList = ampList;                % 冻结扫幅：激励幅值列表（logspace）
cfg.nSeed   = nSeed;                  % 每个幅值下的随机种子个数（统计用）
cfg.doSweepId = doSweepId;            % true=扫频识别共振频率；false=用 DF/理论快路径
cfg.nFreqId = 11;                     % 扫频识别时的频率点数
cfg.T_id    = 20;                     % 扫频每个频点的仿真时长 s
cfg.lambdaList = [0, 0.25, 0.5, 1, 2, 4, 8];  % 自适应场景三：lambda 扫描网格
cfg.tStep = 60;                       % 自适应场景二：幅值阶跃时刻 s
cfg.ampStepGain = 3;                  % 自适应场景二：阶跃后幅值倍数（×3）


steps = options.steps;

%% ===== Step 1: 激励 + 包络诊断 =====
if any(steps == 1)
    fprintf('\n######## Step1: genExcitation + envelopeStats ########\n');
    tCtrl = (0:1/sim.fs_ctrl:sim.T).';
    [xdd, excStats] = genExcitation(exc, tCtrl);

    % 用无控短仿真得到响应包络
    ctlP = ctl; ctlP.adapt_on = false; ctlP.mode = 'passive';
    logP = simDvaOnly(plant, sim, fltSys, ctlP, tCtrl, xdd);
    env = envelopeStats(logP.y, logP.t, sim.t_discard);
    fprintf('  带通时间常数=%.3fs, tau_env=%.3fs\n', tau_bp, env.tau_env);

    fig = figure('Name', 'Step1-Excitation', 'Color', 'w');
    tiledlayout(2,1);
    nexttile;
    plot(tCtrl, xdd); grid on; xlabel('t (s)'); ylabel('xdd');
    title(sprintf('激励 (RMS=%.3e)', excStats.rms));
    nexttile;
    plot(excStats.psd_f, excStats.psd_p); grid on;
    xlim([0, 10]); xlabel('f (Hz)'); ylabel('PSD');
    title('激励功率谱');
    exportgraphics(fig, fullfile(outDir, 'step1_excitation.png'), 'Resolution', 200);
    savefig(fig, fullfile(outDir, 'step1_excitation.fig'));
    save(fullfile(outDir, 'step1_excitation.mat'), 'excStats', 'env', 'xdd', 'tCtrl', 'logP');
end

%% ===== Step 2: 无控幅值-反共振漂移 =====
if any(steps == 2)
    fprintf('\n######## Step2: 无控幅值-频率漂移 ########\n');
    nA = numel(ampList);
    fDf = nan(nA,1);
    Amean = nan(nA,1);
    for i = 1:nA
        excI = exc; excI.amp = ampList(i); excI.seed = 1;
        % 纯简谐更利于复现 backbone
        excI.bw = 0; excI.rand_ratio = 0;
        tCtrl = (0:1/sim.fs_ctrl:sim.T).';
        [xdd, ~] = genExcitation(excI, tCtrl);
        ctlP = ctl; ctlP.adapt_on = false; ctlP.mode = 'passive';
        logP = simDvaOnly(plant, sim, fltSys, ctlP, tCtrl, xdd);
        met = computeMetrics(logP, plant, sim);
        fDf(i) = met.f_df_Hz;
        Amean(i) = met.A_env;
        fprintf('  amp=%.3e, A=%.3e, f_df=%.3f Hz\n', ampList(i), Amean(i), fDf(i));
    end
    fig = figure('Name', 'Step2-Backbone', 'Color', 'w');
    semilogx(ampList, fDf, '-o', 'LineWidth', 1.2); hold on;
    yline(plant.f1_Hz, '--', 'f1 linear');
    grid on; xlabel('amp'); ylabel('f_{DF} (Hz)');
    title('无控: 幅值-等效固有频率（描述函数）');
    exportgraphics(fig, fullfile(outDir, 'step2_backbone.png'), 'Resolution', 200);
    savefig(fig, fullfile(outDir, 'step2_backbone.fig'));
    save(fullfile(outDir, 'step2_backbone.mat'), 'ampList', 'fDf', 'Amean');
end

%% ===== Step 3: 冻结天花板 =====
if any(steps == 3)
    fprintf('\n######## Step3: runFrozenSweep ########\n');
    resultsFrozen = runFrozenSweep(cfg);
    assignin('base', 'resultsFrozen', resultsFrozen);
end

%% ===== Step 4: 自适应三场景 =====
if any(steps == 4)
    fprintf('\n######## Step4: runAdaptive ########\n');
    if options.quick
        cfg.sim.T = max(cfg.sim.T, 80);  % 阶跃需要足够时长
        cfg.tStep = min(40, cfg.sim.T*0.5);
    end
    resultsAdapt = runAdaptive(cfg);
    assignin('base', 'resultsAdapt', resultsAdapt);
end

fprintf('\n[main] 完成。结果目录: %s\n', outDir);
end
