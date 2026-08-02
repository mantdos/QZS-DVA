function results = runAdaptive(cfg)
% runAdaptive  自适应控制器最小验证路径
%
%   results = runAdaptive(cfg)
%
%   当前只跑场景一的 adaptive；passive/frozen/理论 dk/场景二 暂注释。

arguments
    cfg struct
end

plant  = cfg.plant;
host   = cfg.host;
simCfg = cfg.sim;
fltSys = cfg.fltSys;
flt    = cfg.flt;
exc0   = cfg.exc;
ctl0   = cfg.ctl;
outDir = cfg.outDir;

if ~exist(outDir, 'dir'), mkdir(outDir); end

results = struct();
results.scene1 = localScene1(plant, host, simCfg, fltSys, flt, ctl0, exc0, outDir, cfg);
% results.scene2 = localScene2(plant, host, simCfg, fltSys, flt, ctl0, exc0, outDir, cfg);  % 暂不跑阶跃

save(fullfile(outDir, 'runAdaptive_results.mat'), 'results', 'cfg');
fprintf('\n[runAdaptive] 结果已保存至 %s\n', outDir);
end

%% ===== 场景一：固定幅值，只跑 adaptive =====
function S = localScene1(plant, host, simCfg, fltSys, flt, ctl0, exc0, outDir, cfg)
fprintf('\n========== Scene1: 固定幅值（仅 adaptive）==========\n');
tCtrl = (0:1/simCfg.fs_ctrl:simCfg.T).';
[xdd, excStats] = genExcitation(exc0, tCtrl, host, flt);
S.excStats = excStats;

% --- 以下对照 / 理论 dk 暂注释（走通后再打开）---
% ctlP = ctl0; ctlP.adapt_on = false; ctlP.mode = 'passive';
% logP = simDvaOnly(plant, simCfg, fltSys, ctlP, tCtrl, xdd);
% idx = logP.t >= simCfg.t_discard;
% A = mean(abs(hilbert(logP.y(idx))));
% eq = equivStiffnessTheory(plant.k3, plant.k5, A, std(logP.y(idx)));
%
% ctlF = ctl0; ctlF.adapt_on = false; ctlF.mode = 'frozen';
% ctlF.dk_frozen = eq.dk_harmonic; ctlF.dc_frozen = 0;
% logF = simDvaOnly(plant, simCfg, fltSys, ctlF, tCtrl, xdd);

% --- 核心：自适应 ---
ctlA = ctl0; ctlA.adapt_on = true; ctlA.mode = 'adaptive';
logA = simDvaOnly(plant, simCfg, fltSys, ctlA, tCtrl, xdd);

% --- 指标 / 对照打印暂注释 ---
% metP = computeMetrics(logP, plant, simCfg);
% metF = computeMetrics(logF, plant, simCfg);
% metA = computeMetrics(logA, plant, simCfg);
% S.logP = logP; S.logF = logF; S.eq = eq;
% S.metP = metP; S.metF = metF; S.metA = metA;
% if isfield(cfg, 'y_qzs_half')
%     localPrintWorkPoint(logA, cfg.y_qzs_half, simCfg.t_discard);
% end
% fprintf('  e_rms   : pas=%.3e  frz=%.3e  adp=%.3e\n', metP.e_rms, metF.e_rms, metA.e_rms);
% fprintf('  edd_rms : pas=%.3e  frz=%.3e  adp=%.3e\n', metP.edd_rms, metF.edd_rms, metA.edd_rms);
% fprintf('  u_rms   : pas=%.3e  frz=%.3e  adp=%.3e\n', metP.u_rms, metF.u_rms, metA.u_rms);
% fprintf('  dk_hat 稳态=%.4e  vs  理论谐波=%.4e  (高斯=%.4e)\n', ...
%     metA.dk_mean, eq.dk_harmonic, eq.dk_gaussian);

S.logA = logA;
idx = logA.t >= simCfg.t_discard;
fprintf('  [adaptive] e_rms=%.3e, edd_rms=%.3e, u_rms=%.3e, dk_mean=%.4e\n', ...
    rms(logA.e(idx)), rms(logA.edd(idx)), rms(logA.u(idx)), mean(logA.dk_hat(idx)));

% 左列：参考(红实线) vs 实际(黑虚线)；右列：对应跟踪误差(蓝实线)
pos0 = get(0, 'defaultFigurePosition');
fig = figure('Name', 'Track & Error: y / yd / ydd', 'Color', 'w');
tl = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
t = logA.t;

% --- y ---
nexttile; hold on; grid on;
plot(t, logA.ym, 'r-', 'LineWidth', 1.0);
plot(t, logA.y,  'k--', 'LineWidth', 1.0);
ylabel('y (m)');
legend({'参考', '实际'}, 'Location', 'best');

nexttile; hold on; grid on;
plot(t, logA.e, 'b-', 'LineWidth', 1.0);
ylabel('e (m)');
legend({'误差'}, 'Location', 'best');

% --- yd ---
nexttile; hold on; grid on;
plot(t, logA.ym_d, 'r-', 'LineWidth', 1.0);
plot(t, logA.yd,   'k--', 'LineWidth', 1.0);
ylabel('yd (m/s)');
legend({'参考', '实际'}, 'Location', 'best');

nexttile; hold on; grid on;
plot(t, logA.ed, 'b-', 'LineWidth', 1.0);
ylabel('ed (m/s)');
legend({'误差'}, 'Location', 'best');

% --- ydd ---
nexttile; hold on; grid on;
plot(t, logA.ym_dd, 'r-', 'LineWidth', 1.0);
plot(t, logA.ydd,   'k--', 'LineWidth', 1.0);
ylabel('ydd (m/s^2)');
legend({'参考', '实际'}, 'Location', 'best');

nexttile; hold on; grid on;
plot(t, logA.edd, 'b-', 'LineWidth', 1.0);
ylabel('edd (m/s^2)');
legend({'误差'}, 'Location', 'best');

xlabel(tl, 't (s)');
title(tl, '局部量跟踪与误差：参考模型 vs 实际输出');

% 自适应参数与控制力：dk_hat / dc_hat / u
fig2 = figure('Name', 'Adapt & Control: dk / dc / u', 'Color', 'w');
tl2 = tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on; grid on;
plot(t, logA.dk_hat, 'b-', 'LineWidth', 1.0);
ylabel('\Delta k_{hat}');
legend({'\Delta k_{hat}'}, 'Location', 'best');

nexttile; hold on; grid on;
plot(t, logA.dc_hat, 'b-', 'LineWidth', 1.0);
ylabel('\Delta c_{hat}');
legend({'\Delta c_{hat}'}, 'Location', 'best');

nexttile; hold on; grid on;
plot(t, logA.u, 'b-', 'LineWidth', 1.0);
ylabel('u (N)');
legend({'u'}, 'Location', 'best');

xlabel(tl2, 't (s)');
title(tl2, '自适应估计与控制力');

end

%% ===== 场景二：幅值阶跃（暂注释，整函数保留）=====
function S = localScene2(plant, host, simCfg, fltSys, flt, ctl0, exc0, outDir, cfg)
% 当前未从入口调用；需要时取消上面 results.scene2 = ... 的注释即可。
fprintf('\n========== Scene2: 幅值阶跃（当前未启用）==========\n');
S = struct();
% tStep = 30;
% if isfield(cfg, 'tStep'), tStep = cfg.tStep; end
% gain = 3;
% if isfield(cfg, 'ampStepGain'), gain = cfg.ampStepGain; end
%
% tCtrl = (0:1/simCfg.fs_ctrl:simCfg.T).';
% [xdd0, ~] = genExcitation(exc0, tCtrl, host, flt);
% exc2 = exc0; exc2.amp = exc0.amp * gain;
% [xdd1, ~] = genExcitation(exc2, tCtrl, host, flt);
% xdd = xdd0;
% xdd(tCtrl >= tStep) = xdd1(tCtrl >= tStep);
%
% ctlP0 = ctl0; ctlP0.adapt_on = false; ctlP0.mode = 'passive';
% logP0 = simDvaOnly(plant, simCfg, fltSys, ctlP0, tCtrl, xdd0);
% idx0 = logP0.t >= simCfg.t_discard;
% A0 = mean(abs(hilbert(logP0.y(idx0))));
% eq0 = equivStiffnessTheory(plant.k3, plant.k5, A0, std(logP0.y(idx0)));
%
% ctlPas = ctl0; ctlPas.adapt_on = false; ctlPas.mode = 'passive';
% ctlFrz = ctl0; ctlFrz.adapt_on = false; ctlFrz.mode = 'frozen';
% ctlFrz.dk_frozen = eq0.dk_harmonic; ctlFrz.dc_frozen = 0;
% ctlAdp = ctl0; ctlAdp.adapt_on = true; ctlAdp.mode = 'adaptive';
%
% logPas = simDvaOnly(plant, simCfg, fltSys, ctlPas, tCtrl, xdd);
% logFrz = simDvaOnly(plant, simCfg, fltSys, ctlFrz, tCtrl, xdd);
% logAdp = simDvaOnly(plant, simCfg, fltSys, ctlAdp, tCtrl, xdd);
% ... 其余绘图 / 重收敛统计见 git 历史
end

%% ------------------------------------------------------------------------
% function tRec = localReconvTime(t, dk, tStep, dkFinal)
% % 进入终值 ±10% 带并保持到结束的时刻（相对阶跃时刻）
% band = 0.1 * abs(dkFinal);
% if band < eps
%     band = 0.1 * max(abs(dk(t >= tStep)));
% end
% if ~isfinite(band) || band < eps
%     tRec = Inf;
%     return;
% end
% idx = t >= tStep;
% tt = t(idx);
% yy = dk(idx);
% inBand = abs(yy - dkFinal) <= band;
% lastOut = find(~inBand, 1, 'last');
% if isempty(lastOut)
%     tRec = 0;
% elseif lastOut == numel(inBand)
%     tRec = Inf;
% else
%     tRec = tt(lastOut + 1) - tStep;
% end
% end
%
% function localPrintWorkPoint(logP, y_qzs_half, tDiscard)
% idx = logP.t >= tDiscard;
% yPeak = max(abs(logP.y(idx)));
% ratio = yPeak / y_qzs_half;
% fprintf('  QZS 工作点: max|y|/y_qzs_half = %.3f', ratio);
% if ratio >= 1
%     fprintf('  ← 已离开低刚度区，请降低 exc.amp');
% end
% fprintf('\n');
% end
