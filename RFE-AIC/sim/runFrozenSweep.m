function results = runFrozenSweep(cfg)
% runFrozenSweep  参数冻结天花板测试（第一优先实验）
%
%   results = runFrozenSweep(cfg)
%
%   流程:
%     1) 扫描激励幅值 amp（logspace）
%     2) 无控仿真 → 统计 σ_y、A → 理论 dk (谐波/高斯)
%     3) 冻结参数受控仿真（版本 A / B）
%     4) 多随机种子取统计；每幅值做一次扫频识别共振频率
%     5) 出图并保存 .mat / .fig / .png

arguments
    cfg struct
end

plant  = cfg.plant;
simCfg = cfg.sim;
fltSys = cfg.fltSys;
exc0   = cfg.exc;
ctl0   = cfg.ctl;
outDir = cfg.outDir;

ampList = cfg.ampList(:);
nSeed   = cfg.nSeed;
nAmp    = numel(ampList);
doSweepId = ~isfield(cfg, 'doSweepId') || cfg.doSweepId;
nFreqId = 11;
if isfield(cfg, 'nFreqId'), nFreqId = cfg.nFreqId; end
T_id = 20;
if isfield(cfg, 'T_id'), T_id = cfg.T_id; end

R.amp = ampList;
R.f_passive   = nan(nAmp, 1);
R.f_frozenA   = nan(nAmp, 1);
R.f_frozenB   = nan(nAmp, 1);
R.f_df_pas    = nan(nAmp, nSeed);
R.yrms_pas    = nan(nAmp, nSeed);
R.yrms_frA    = nan(nAmp, nSeed);
R.yrms_frB    = nan(nAmp, nSeed);
R.A_cv_pas    = nan(nAmp, nSeed);
R.dkA         = nan(nAmp, nSeed);
R.dkB         = nan(nAmp, nSeed);
R.e_rms_frA   = nan(nAmp, nSeed);
R.edd_rms_frA = nan(nAmp, nSeed);
R.e_rms_frB   = nan(nAmp, nSeed);
R.edd_rms_frB = nan(nAmp, nSeed);
R.u_rms_frA   = nan(nAmp, nSeed);

fprintf('\n========== runFrozenSweep: %d amps × %d seeds ==========\n', nAmp, nSeed);

for iA = 1:nAmp
    amp = ampList(iA);
    fprintf('\n--- amp = %.4e (%d/%d) ---\n', amp, iA, nAmp);

    for iS = 1:nSeed
        fprintf('  seed %d/%d\n', iS, nSeed);
        exc = exc0;
        exc.amp  = amp;
        exc.seed = iS;

        tCtrl = (0:1/simCfg.fs_ctrl:simCfg.T).';
        [xdd, ~] = genExcitation(exc, tCtrl);

        % --- 无控 ---
        ctlP = ctl0;
        ctlP.adapt_on = false;
        ctlP.mode = 'passive';
        logP = simDvaOnly(plant, simCfg, fltSys, ctlP, tCtrl, xdd);
        envP = envelopeStats(logP.y, logP.t, simCfg.t_discard);
        metP = computeMetrics(logP, plant, simCfg);

        idx = logP.t >= simCfg.t_discard;
        eq = equivStiffnessTheory(plant.k3, plant.k5, envP.A_mean, std(logP.y(idx)));

        R.dkA(iA,iS)      = eq.dk_harmonic;
        R.dkB(iA,iS)      = eq.dk_gaussian;
        R.yrms_pas(iA,iS) = metP.y_rms;
        R.A_cv_pas(iA,iS) = envP.A_cv;
        R.f_df_pas(iA,iS) = metP.f_df_Hz;

        % --- 冻结 A ---
        ctlA = ctl0;
        ctlA.adapt_on  = false;
        ctlA.mode      = 'frozen';
        ctlA.dk_frozen = eq.dk_harmonic;
        ctlA.dc_frozen = 0;
        logA = simDvaOnly(plant, simCfg, fltSys, ctlA, tCtrl, xdd);
        metA = computeMetrics(logA, plant, simCfg);
        R.yrms_frA(iA,iS)    = metA.y_rms;
        R.e_rms_frA(iA,iS)   = metA.e_rms;
        R.edd_rms_frA(iA,iS) = metA.edd_rms;
        R.u_rms_frA(iA,iS)   = metA.u_rms;

        % --- 冻结 B ---
        ctlB = ctl0;
        ctlB.adapt_on  = false;
        ctlB.mode      = 'frozen';
        ctlB.dk_frozen = eq.dk_gaussian;
        ctlB.dc_frozen = 0;
        logB = simDvaOnly(plant, simCfg, fltSys, ctlB, tCtrl, xdd);
        metB = computeMetrics(logB, plant, simCfg);
        R.yrms_frB(iA,iS)    = metB.y_rms;
        R.e_rms_frB(iA,iS)   = metB.e_rms;
        R.edd_rms_frB(iA,iS) = metB.edd_rms;
    end

    % --- 每幅值一次扫频识别（用多种子平均 dk）---
    dkA_mu = mean(R.dkA(iA,:), 'omitnan');
    dkB_mu = mean(R.dkB(iA,:), 'omitnan');
    fprintf('  dkA=%.4e, dkB=%.4e\n', dkA_mu, dkB_mu);

    ctlP = ctl0; ctlP.adapt_on = false; ctlP.mode = 'passive';
    ctlA = ctl0; ctlA.adapt_on = false; ctlA.mode = 'frozen';
    ctlA.dk_frozen = dkA_mu; ctlA.dc_frozen = 0;
    ctlB = ctl0; ctlB.adapt_on = false; ctlB.mode = 'frozen';
    ctlB.dk_frozen = dkB_mu; ctlB.dc_frozen = 0;

    if doSweepId
        fprintf('  扫频识别共振频率...\n');
        R.f_passive(iA) = identifyResonanceFreq(plant, simCfg, fltSys, ctlP, amp, ...
            'nFreq', nFreqId, 'T', T_id);
        R.f_frozenA(iA) = identifyResonanceFreq(plant, simCfg, fltSys, ctlA, amp, ...
            'nFreq', nFreqId, 'T', T_id);
        R.f_frozenB(iA) = identifyResonanceFreq(plant, simCfg, fltSys, ctlB, amp, ...
            'nFreq', nFreqId, 'T', T_id);
    else
        % 快速路径：无控用 DF；受控报线性固有频率（理想天花板）
        R.f_passive(iA) = mean(R.f_df_pas(iA,:), 'omitnan');
        R.f_frozenA(iA) = plant.f1_Hz;
        R.f_frozenB(iA) = plant.f1_Hz;
    end
    fprintf('  f: pas=%.3f, frA=%.3f, frB=%.3f Hz\n', ...
        R.f_passive(iA), R.f_frozenA(iA), R.f_frozenB(iA));
end

R.f1_lin_Hz = plant.f1_Hz;
results = R;

localPlotFrozen(R, outDir);
save(fullfile(outDir, 'runFrozenSweep_results.mat'), 'results', 'cfg');
fprintf('\n[runFrozenSweep] 结果已保存至 %s\n', outDir);
end

%% ------------------------------------------------------------------------
function localPlotFrozen(R, outDir)
amp = R.amp;

fig1 = figure('Name', 'FrozenSweep-Fig1', 'Color', 'w');
tiledlayout(1,1); nexttile;
semilogx(amp, R.f_passive, '-o', 'LineWidth', 1.2); hold on;
semilogx(amp, R.f_frozenA, '-s', 'LineWidth', 1.2);
semilogx(amp, R.f_frozenB, '-d', 'LineWidth', 1.2);
yline(R.f1_lin_Hz, '--', '线性固有频率');
grid on;
xlabel('激励幅值 amp'); ylabel('反共振/共振频率 (Hz)');
legend({'无控','冻结A(谐波)','冻结B(高斯)'}, 'Location', 'best');
title('天花板测试: 幅值-反共振频率');
exportgraphics(fig1, fullfile(outDir, 'frozen_fig1_f_ar.png'), 'Resolution', 200);
savefig(fig1, fullfile(outDir, 'frozen_fig1_f_ar.fig'));

fig2 = figure('Name', 'FrozenSweep-Fig2', 'Color', 'w');
tiledlayout(1,1); nexttile;
loglog(amp, mean(R.yrms_pas,2,'omitnan'), '-o', 'LineWidth', 1.2); hold on;
loglog(amp, mean(R.yrms_frA,2,'omitnan'), '-s', 'LineWidth', 1.2);
loglog(amp, mean(R.yrms_frB,2,'omitnan'), '-d', 'LineWidth', 1.2);
grid on;
xlabel('激励幅值 amp'); ylabel('y RMS (m)');
legend({'无控','冻结A','冻结B'}, 'Location', 'best');
title('响应 RMS vs 激励幅值');
exportgraphics(fig2, fullfile(outDir, 'frozen_fig2_yrms.png'), 'Resolution', 200);
savefig(fig2, fullfile(outDir, 'frozen_fig2_yrms.fig'));

fig3 = figure('Name', 'FrozenSweep-Fig3', 'Color', 'w');
tiledlayout(1,1); nexttile;
errorbar(amp, mean(R.A_cv_pas,2,'omitnan')*100, ...
    std(R.A_cv_pas,0,2,'omitnan')*100, '-o', 'LineWidth', 1.2);
set(gca, 'XScale', 'log'); grid on;
xlabel('激励幅值 amp'); ylabel('A_{cv} (%)');
yline(10, '--'); yline(30, '--');
title('包络变异系数 A_{cv}');
exportgraphics(fig3, fullfile(outDir, 'frozen_fig3_Acv.png'), 'Resolution', 200);
savefig(fig3, fullfile(outDir, 'frozen_fig3_Acv.fig'));
end
