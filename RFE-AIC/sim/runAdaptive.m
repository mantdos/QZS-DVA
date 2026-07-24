function results = runAdaptive(cfg)
% runAdaptive  自适应性能测试（天花板确认后运行）
%
%   results = runAdaptive(cfg)
%
%   场景一: 固定幅值，观察参数收敛 / 漂移
%   场景二: 幅值阶跃（默认 t=60s 时 ×3）
%   场景三: lambda 扫描，检查中间最优

arguments
    cfg struct
end

plant  = cfg.plant;
simCfg = cfg.sim;
fltSys = cfg.fltSys;
exc0   = cfg.exc;
ctl0   = cfg.ctl;
outDir = cfg.outDir;

if ~exist(outDir, 'dir'), mkdir(outDir); end

results = struct();
results.scene1 = localScene1(plant, simCfg, fltSys, ctl0, exc0, outDir);
results.scene2 = localScene2(plant, simCfg, fltSys, ctl0, exc0, outDir, cfg);
results.scene3 = localScene3(plant, simCfg, fltSys, ctl0, exc0, outDir, cfg);

save(fullfile(outDir, 'runAdaptive_results.mat'), 'results', 'cfg');
fprintf('\n[runAdaptive] 结果已保存至 %s\n', outDir);
end

%% ===== 场景一：固定幅值收敛 =====
function S = localScene1(plant, simCfg, fltSys, ctl0, exc0, outDir)
fprintf('\n========== Scene1: 固定幅值自适应 ==========\n');
exc = exc0;
tCtrl = (0:1/simCfg.fs_ctrl:simCfg.T).';
[xdd, ~] = genExcitation(exc, tCtrl);

% 无控 → 理论 dk 参考
ctlP = ctl0; ctlP.adapt_on = false; ctlP.mode = 'passive';
logP = simDvaOnly(plant, simCfg, fltSys, ctlP, tCtrl, xdd);
envP = envelopeStats(logP.y, logP.t, simCfg.t_discard);
idx = logP.t >= simCfg.t_discard;
eq = equivStiffnessTheory(plant.k3, plant.k5, envP.A_mean, std(logP.y(idx)));

% 冻结最优（谐波）对照
ctlF = ctl0; ctlF.adapt_on = false; ctlF.mode = 'frozen';
ctlF.dk_frozen = eq.dk_harmonic; ctlF.dc_frozen = 0;
logF = simDvaOnly(plant, simCfg, fltSys, ctlF, tCtrl, xdd);

% 自适应
ctlA = ctl0; ctlA.adapt_on = true; ctlA.mode = 'adaptive';
logA = simDvaOnly(plant, simCfg, fltSys, ctlA, tCtrl, xdd);
envA = envelopeStats(logA.y, logA.t, simCfg.t_discard);
metA = computeMetrics(logA, plant, simCfg);
metF = computeMetrics(logF, plant, simCfg);

S.logA = logA; S.logF = logF; S.eq = eq; S.metA = metA; S.metF = metF; S.envA = envA;

fig = figure('Name', 'Adaptive-Scene1', 'Color', 'w');
tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(logA.t, logA.dk_hat, 'LineWidth', 1); hold on;
yline(eq.dk_harmonic, '--r', 'dk_A');
yline(eq.dk_gaussian, '--', 'dk_B');
grid on; xlabel('t (s)'); ylabel('\Delta k_{hat}');
title('刚度自适应（观察漂移）'); legend({'\Delta k_{hat}','谐波理论','高斯理论'});

nexttile;
plot(logA.t, logA.dc_hat, 'LineWidth', 1); grid on;
xlabel('t (s)'); ylabel('\Delta c_{hat}'); title('阻尼自适应');

nexttile;
plot(logA.t, logA.s, 'LineWidth', 0.8); grid on;
xlabel('t (s)'); ylabel('s'); title('复合误差面');

nexttile;
plot(logA.t, logA.e, 'LineWidth', 0.8); grid on;
xlabel('t (s)'); ylabel('e'); title(sprintf('跟踪误差 (e_{rms}=%.2e)', metA.e_rms));

title(tl, '场景一: 固定幅值自适应');
exportgraphics(fig, fullfile(outDir, 'adapt_scene1.png'), 'Resolution', 200);
savefig(fig, fullfile(outDir, 'adapt_scene1.fig'));

fprintf('  dk_A=%.4e, dk_B=%.4e, dk_hat(稳态)=%.4e\n', ...
    eq.dk_harmonic, eq.dk_gaussian, metA.dk_mean);
fprintf('  e_rms: adapt=%.3e, frozen=%.3e | edd_rms: adapt=%.3e, frozen=%.3e\n', ...
    metA.e_rms, metF.e_rms, metA.edd_rms, metF.edd_rms);
end

%% ===== 场景二：幅值阶跃 =====
function S = localScene2(plant, simCfg, fltSys, ctl0, exc0, outDir, cfg)
fprintf('\n========== Scene2: 幅值阶跃 ==========\n');
tStep = 60;
if isfield(cfg, 'tStep'), tStep = cfg.tStep; end
gain = 3;
if isfield(cfg, 'ampStepGain'), gain = cfg.ampStepGain; end

exc = exc0;
tCtrl = (0:1/simCfg.fs_ctrl:simCfg.T).';
[xdd0, ~] = genExcitation(exc, tCtrl);
amp2 = exc.amp * gain;
exc2 = exc; exc2.amp = amp2;
[xdd1, ~] = genExcitation(exc2, tCtrl);
xdd = xdd0;
xdd(tCtrl >= tStep) = xdd1(tCtrl >= tStep);

% 初始幅值下的冻结 dk
ctlP = ctl0; ctlP.adapt_on = false; ctlP.mode = 'passive';
logP0 = simDvaOnly(plant, simCfg, fltSys, ctlP, tCtrl, xdd0);
env0 = envelopeStats(logP0.y, logP0.t, simCfg.t_discard);
idx = logP0.t >= simCfg.t_discard;
eq0 = equivStiffnessTheory(plant.k3, plant.k5, env0.A_mean, std(logP0.y(idx)));

% 三种模式
ctlPas = ctl0; ctlPas.adapt_on = false; ctlPas.mode = 'passive';
ctlFrz = ctl0; ctlFrz.adapt_on = false; ctlFrz.mode = 'frozen';
ctlFrz.dk_frozen = eq0.dk_harmonic; ctlFrz.dc_frozen = 0;
ctlAdp = ctl0; ctlAdp.adapt_on = true; ctlAdp.mode = 'adaptive';

logPas = simDvaOnly(plant, simCfg, fltSys, ctlPas, tCtrl, xdd);
logFrz = simDvaOnly(plant, simCfg, fltSys, ctlFrz, tCtrl, xdd);
logAdp = simDvaOnly(plant, simCfg, fltSys, ctlAdp, tCtrl, xdd);

S.logPas = logPas; S.logFrz = logFrz; S.logAdp = logAdp;
S.eq0 = eq0; S.tStep = tStep; S.gain = gain;

% 阶跃后指标
idx2 = tCtrl >= tStep + simCfg.t_discard;
S.yrms_pas = rms(logPas.y(idx2));
S.yrms_frz = rms(logFrz.y(idx2));
S.yrms_adp = rms(logAdp.y(idx2));
S.e_rms_adp = rms(logAdp.e(idx2));
S.edd_rms_adp = rms(logAdp.edd(idx2));

fig = figure('Name', 'Adaptive-Scene2', 'Color', 'w');
tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
plot(tCtrl, logPas.y, tCtrl, logFrz.y, tCtrl, logAdp.y, 'LineWidth', 0.7);
xline(tStep, '--k'); grid on;
xlabel('t (s)'); ylabel('y'); legend({'无控','冻结','自适应'});
title('位移响应');

nexttile;
plot(tCtrl, logAdp.dk_hat, 'LineWidth', 1); hold on;
yline(eq0.dk_harmonic, '--r', '初始冻结值');
xline(tStep, '--k'); grid on;
xlabel('t (s)'); ylabel('\Delta k_{hat}'); title('自适应刚度');

nexttile;
plot(tCtrl, logAdp.e, 'LineWidth', 0.7); xline(tStep, '--k'); grid on;
xlabel('t (s)'); ylabel('e'); title('自适应跟踪误差');

nexttile;
plot(tCtrl, abs(hilbert(logPas.y)), tCtrl, abs(hilbert(logAdp.y)), 'LineWidth', 0.8);
xline(tStep, '--k'); grid on;
xlabel('t (s)'); ylabel('|env(y)|'); legend({'无控','自适应'});
title(sprintf('包络 (阶跃后 y_{rms}: pas=%.2e adp=%.2e)', S.yrms_pas, S.yrms_adp));

title(tl, sprintf('场景二: t=%.0fs 幅值 ×%.1f', tStep, gain));
exportgraphics(fig, fullfile(outDir, 'adapt_scene2.png'), 'Resolution', 200);
savefig(fig, fullfile(outDir, 'adapt_scene2.fig'));
end

%% ===== 场景三：lambda 扫描 =====
function S = localScene3(plant, simCfg, fltSys, ctl0, exc0, outDir, cfg)
fprintf('\n========== Scene3: lambda 扫描 ==========\n');
if isfield(cfg, 'lambdaList')
    lamList = cfg.lambdaList(:);
else
    lamList = [0, 0.25, 0.5, 1, 2, 4, 8];
end
nL = numel(lamList);
S.lambda = lamList;
S.e_rms   = nan(nL,1);
S.edd_rms = nan(nL,1);
S.y_rms   = nan(nL,1);
S.u_rms   = nan(nL,1);
S.u_over_Fpas = nan(nL,1);
S.dk_mean = nan(nL,1);

exc = exc0;
tCtrl = (0:1/simCfg.fs_ctrl:simCfg.T).';
[xdd, ~] = genExcitation(exc, tCtrl);
wc = fltSys.omega_c;

for i = 1:nL
    lam = lamList(i);
    fprintf('  lambda = %.3f (%d/%d)\n', lam, i, nL);
    ctl = ctl0;
    ctl.adapt_on = true;
    ctl.mode = 'adaptive';
    ctl.lambda = lam;
    ctl.c2 = lam / wc;
    ctl.Ms = 1 + lam;

    logA = simDvaOnly(plant, simCfg, fltSys, ctl, tCtrl, xdd);
    met = computeMetrics(logA, plant, simCfg);
    S.e_rms(i) = met.e_rms;
    S.edd_rms(i) = met.edd_rms;
    S.y_rms(i) = met.y_rms;
    S.u_rms(i) = met.u_rms;
    S.u_over_Fpas(i) = met.u_over_Fpas;
    S.dk_mean(i) = met.dk_mean;
    envelopeStats(logA.y, logA.t, simCfg.t_discard);
end

fig = figure('Name', 'Adaptive-Scene3', 'Color', 'w');
tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
plot(lamList, S.e_rms, '-o', 'LineWidth', 1.2); grid on;
xlabel('\lambda'); ylabel('||e||_{rms}'); title('状态误差');
nexttile;
plot(lamList, S.edd_rms, '-o', 'LineWidth', 1.2); grid on;
xlabel('\lambda'); ylabel('||?||_{rms}'); title('加速度误差（检验反作用力项）');
nexttile;
plot(lamList, S.y_rms, '-o', 'LineWidth', 1.2); grid on;
xlabel('\lambda'); ylabel('y_{rms}'); title('位移 RMS');
nexttile;
plot(lamList, S.u_over_Fpas, '-o', 'LineWidth', 1.2); grid on;
xlabel('\lambda'); ylabel('u_{rms}/F_{pas,rms}'); title('控制力相对被动恢复力');
title(tl, '场景三: \lambda 扫描（是否存在中间最优）');
exportgraphics(fig, fullfile(outDir, 'adapt_scene3_lambda.png'), 'Resolution', 200);
savefig(fig, fullfile(outDir, 'adapt_scene3_lambda.fig'));
end
