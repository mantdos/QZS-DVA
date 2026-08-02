function results = runLambdaSweep(cfg)
% runLambdaSweep  lambda 扫描（独立脚本，本轮不在 main 中调用）
%
%   results = runLambdaSweep(cfg)
%
%   用途: 后续验证 λ 价值时单独运行。
%   要求 cfg 含 plant/host/sim/fltSys/flt/ctl/exc/outDir，
%   可选 cfg.lambdaList。

arguments
    cfg struct
end

plant  = cfg.plant;
host   = cfg.host;
simCfg = cfg.sim;
fltSys = cfg.fltSys;
flt    = cfg.flt;
ctl0   = cfg.ctl;
exc0   = cfg.exc;
outDir = cfg.outDir;

if isfield(cfg, 'lambdaList')
    lamList = cfg.lambdaList(:);
else
    lamList = [0, 0.25, 0.5, 1, 2, 4, 8];
end

if ~exist(outDir, 'dir'), mkdir(outDir); end

nL = numel(lamList);
results.lambda      = lamList;
results.e_rms       = nan(nL, 1);
results.edd_rms     = nan(nL, 1);
results.y_rms       = nan(nL, 1);
results.u_rms       = nan(nL, 1);
results.u_over_Fpas = nan(nL, 1);
results.dk_mean     = nan(nL, 1);

tCtrl = (0:1/simCfg.fs_ctrl:simCfg.T).';
[xdd, ~] = genExcitation(exc0, tCtrl, host, flt);
wc = fltSys.omega_c;

fprintf('\n========== runLambdaSweep ==========\n');
for i = 1:nL
    lam = lamList(i);
    fprintf('  lambda = %.3f (%d/%d)\n', lam, i, nL);
    ctl = ctl0;
    ctl.adapt_on = true;
    ctl.mode = 'adaptive';
    ctl.lambda = lam;
    ctl.w2 = lam / wc;
    ctl.Ms = 1 + lam;

    logA = simDvaOnly(plant, simCfg, fltSys, ctl, tCtrl, xdd);
    met = computeMetrics(logA, plant, simCfg);
    results.e_rms(i) = met.e_rms;
    results.edd_rms(i) = met.edd_rms;
    results.y_rms(i) = met.y_rms;
    results.u_rms(i) = met.u_rms;
    results.u_over_Fpas(i) = met.u_over_Fpas;
    results.dk_mean(i) = met.dk_mean;
end

fig = figure('Name', 'LambdaSweep', 'Color', 'w');
tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
plot(lamList, results.e_rms, '-o', 'LineWidth', 1.2); grid on;
xlabel('\lambda'); ylabel('||e||_{rms}'); title('状态误差');
nexttile;
plot(lamList, results.edd_rms, '-o', 'LineWidth', 1.2); grid on;
xlabel('\lambda'); ylabel('||e_{dd}||_{rms}'); title('加速度误差');
nexttile;
plot(lamList, results.y_rms, '-o', 'LineWidth', 1.2); grid on;
xlabel('\lambda'); ylabel('y_{rms}'); title('位移 RMS');
nexttile;
plot(lamList, results.u_over_Fpas, '-o', 'LineWidth', 1.2); grid on;
xlabel('\lambda'); ylabel('u_{rms}/F_{pas,rms}'); title('控制力相对被动恢复力');
title(tl, '\lambda 扫描');
exportgraphics(fig, fullfile(outDir, 'adapt_lambda_sweep.png'), 'Resolution', 200);

save(fullfile(outDir, 'runLambdaSweep_results.mat'), 'results', 'cfg');
fprintf('[runLambdaSweep] 结果已保存至 %s\n', outDir);
end
