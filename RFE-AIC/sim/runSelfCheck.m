function results = runSelfCheck(cfg)
% runSelfCheck  耦合宿主 + 带限激励 物理自检
%
%   results = runSelfCheck(cfg)
%
%   检查项:
%     1) genExcitation: 带限 PSD、总 RMS≈amp、DC 可忽略
%     2) 改变 fs_ctrl 时带内功率基本不变（fs 无关性）
%     3) 线性无控耦合: 宿主传递率双峰+中间反共振 ≈ f_m
%     4) mu→0 时 simCoupled ≈ simHostOnly
%     5) QZS 工作点 max|y| / y_qzs_half

arguments
    cfg struct
end

plant = cfg.plant;
host  = cfg.host;
flt   = cfg.flt;
exc0  = cfg.exc;
sim0  = cfg.sim;

% 自检用较短时长
sim = sim0;
sim.T = min(sim0.T, 40);
sim.t_discard = min(sim0.t_discard, 8);
if isfield(cfg, 'y_qzs_half')
    y_qzs = cfg.y_qzs_half;
else
    y_qzs = sqrt(4*plant.k1 / (3*plant.k3));
end

results = struct();
fprintf('\n========== runSelfCheck ==========\n');

%% 1) 激励带限与 RMS
fprintf('\n--- Check1: genExcitation 带限 / RMS ---\n');
exc = exc0;
exc.mode = 'force';
exc.rand_ratio = max(exc.rand_ratio, 0.15);
exc.amp = 1.0;
tCtrl = (0:1/sim.fs_ctrl:sim.T).';
[F, st] = genExcitation(exc, tCtrl, flt);
results.check1.rms = st.rms;
results.check1.r_inband = st.r_inband;
rmsErr = abs(st.rms - exc.amp) / exc.amp;
fprintf('  RMS=%.4e vs amp=%.4e, 相对误差=%.2f%%\n', st.rms, exc.amp, 100*rmsErr);

% DC / 低频功率占比、通带功率占比
f = st.psd_f(:); p = st.psd_p(:);
Ptot = trapz(f, p);
maskDc = f <= exc.f_lo*0.5;
Pdc = trapz(f(maskDc), p(maskDc));
fprintf('  P(f<%.2f Hz)/Ptot = %.3e （应可忽略）\n', exc.f_lo*0.5, Pdc/Ptot);
maskIn = f >= exc.f_lo & f <= exc.f_hi;
Pin = trapz(f(maskIn), p(maskIn));
fprintf('  通带功率占比 Pin/Ptot = %.3f\n', Pin/Ptot);
results.check1.rms_err = rmsErr;
results.check1.Pdc_frac = Pdc/Ptot;
results.check1.Pin_frac = Pin/Ptot;

%% 2) fs 无关性（带内随机功率占比）
% 不同 fs 下 randn 长度不同 → 实现不同；比较 Pinband/Ptotal 更合理
fprintf('\n--- Check2: fs_ctrl 2000→4000 带内功率占比 ---\n');
exc2 = exc; exc2.seed = 7; exc2.rand_ratio = 1.0;  % 纯随机
t2k = (0:1/2000:sim.T).';
t4k = (0:1/4000:sim.T).';
[~, st2] = genExcitation(exc2, t2k, flt);
[~, st4] = genExcitation(exc2, t4k, flt);
frac2 = localInbandFrac(st2, exc2.f_lo, exc2.f_hi, flt);
frac4 = localInbandFrac(st4, exc2.f_lo, exc2.f_hi, flt);
results.check2.frac_2k = frac2;
results.check2.frac_4k = frac4;
rel = abs(frac4 - frac2) / frac2;
fprintf('  带内功率占比 @2k=%.4f, @4k=%.4f, 相对差=%.2f%%\n', frac2, frac4, 100*rel);
fprintf('  （理想平坦带通 ≈ Δf/(f_hi-f_lo) = %.4f）\n', ...
    (flt.zeta_b*flt.omega_d/pi) / (exc2.f_hi - exc2.f_lo));
results.check2.rel_diff = rel;

%% 3) 线性无控耦合：双峰 + 反共振 ≈ f_m
fprintf('\n--- Check3: 线性无控耦合传递率（反共振） ---\n');
plantL = plant;
plantL.k3 = 0; plantL.k5 = 0;
ctlP = cfg.ctl; ctlP.adapt_on = false; ctlP.mode = 'passive';
fltSys = cfg.fltSys;

% 简谐扫频估计 |X|/|F|，在 f1 邻域找局部最小（非全局最小：高频 |H|~1/ω?）
fGrid = linspace(0.5, 6.0, 41);
TX = zeros(size(fGrid));
simSweep = sim; simSweep.T = 25; simSweep.t_discard = 10;
ampF = 1.0;
for i = 1:numel(fGrid)
    tt = (0:1/simSweep.fs_ctrl:simSweep.T).';
    Fi = ampF * sin(2*pi*fGrid(i)*tt);
    logC = simCoupled(plantL, host, simSweep, fltSys, ctlP, tt, Fi);
    idx = logC.t >= simSweep.t_discard;
    TX(i) = rms(logC.x(idx)) / rms(Fi(idx));
end
% 在 [0.7*f1, 1.3*f1] 内找局部最小
f1 = plant.f1_Hz;
mask = fGrid >= 0.7*f1 & fGrid <= 1.3*f1;
fLoc = fGrid(mask); TLoc = TX(mask);
[~, iMinLoc] = min(TLoc);
iMin = find(mask, 1, 'first') + iMinLoc - 1;
fAr = fGrid(iMin);
f0 = fGrid(iMin-1); f1p = fGrid(iMin); f2 = fGrid(iMin+1);
y0 = TX(iMin-1); y1 = TX(iMin); y2 = TX(iMin+1);
fAr = f1p + 0.5*(y0 - y2)/(y0 - 2*y1 + y2) * (f2 - f0)/2;
% 粗检双峰：左侧与右侧是否存在高于反共振的峰
hasLeftPeak  = any(TX(fGrid < fAr) > TX(iMin) * 2);
hasRightPeak = any(TX(fGrid > fAr) > TX(iMin) * 2);
results.check3.f_antires = fAr;
results.check3.f1 = f1;
results.check3.TX = TX;
results.check3.fGrid = fGrid;
results.check3.has_double_peak = hasLeftPeak && hasRightPeak;
fprintf('  反共振 ≈ %.4f Hz, f_m(f1) = %.4f Hz, |差|=%.4f Hz\n', ...
    fAr, f1, abs(fAr - f1));
fprintf('  双峰结构: 左峰=%d, 右峰=%d\n', hasLeftPeak, hasRightPeak);

%% 4) mu→0：耦合 ≈ 宿主 alone
% 同步缩小 m,c,k*（保持 ω 不变），使 F_react→0 且数值仍稳定
fprintf('\n--- Check4: mu→0 时 coupled ≈ hostOnly ---\n');
epsMu = 1e-4;
plantTiny = plantL;
plantTiny.m   = plantL.m   * epsMu;
plantTiny.c   = plantL.c   * epsMu;
plantTiny.k1  = plantL.k1  * epsMu;
plantTiny.k3  = 0;
plantTiny.k5  = 0;
plantTiny.k_m = plantL.k_m * epsMu;
plantTiny.c_m = plantL.c_m * epsMu;
hostTiny = host;
hostTiny.mu = plantTiny.m / hostTiny.M;
tt = (0:1/sim.fs_ctrl:min(sim.T, 20)).';
excT = exc; excT.rand_ratio = 0; excT.amp = 1.0;
[Fi, ~] = genExcitation(excT, tt, flt);
simT = sim; simT.T = tt(end);
logH = simHostOnly(hostTiny, simT, tt, Fi);
logC = simCoupled(plantTiny, hostTiny, simT, fltSys, ctlP, tt, Fi);
idx = tt >= min(5, sim.t_discard);
relX = rms(logC.x(idx) - logH.x(idx)) / rms(logH.x(idx));
results.check4.rel_x_diff = relX;
results.check4.mu = hostTiny.mu;
fprintf('  mu=%.3e, ||x_c - x_h||_rms / ||x_h||_rms = %.3e\n', hostTiny.mu, relX);

%% 5) QZS 工作点
fprintf('\n--- Check5: QZS 工作点 ---\n');
excW = exc0; excW.mode = 'force';
tt = (0:1/sim.fs_ctrl:sim.T).';
[Fi, ~] = genExcitation(excW, tt, flt);
logW = simCoupled(plant, host, sim, fltSys, ctlP, tt, Fi);
idx = logW.t >= sim.t_discard;
yPeak = max(abs(logW.y(idx)));
ratio = yPeak / y_qzs;
results.check5.y_peak = yPeak;
results.check5.y_qzs = y_qzs;
results.check5.ratio = ratio;
fprintf('  amp=%.3g N(RMS), max|y|=%.3e, y_qzs=%.3e, 比值=%.3f\n', ...
    excW.amp, yPeak, y_qzs, ratio);
if ratio >= 1
    fprintf('  警告: 已离开/接近 QZS 有效区间，请降低 amp\n');
else
    fprintf('  工作点在 QZS 低刚度区内\n');
end

fprintf('\n========== runSelfCheck 完成 ==========\n');
end

%% ------------------------------------------------------------------------
function frac = localInbandFrac(st, f_lo, f_hi, flt)
% 控制器带内功率 / 激励通带总功率
f = st.psd_f(:); p = st.psd_p(:);
dF = flt.zeta_b * flt.omega_d / (2*pi);
f_d = flt.omega_d / (2*pi);
maskTot = f >= f_lo & f <= f_hi;
maskIn  = f >= (f_d - dF) & f <= (f_d + dF);
Ptot = trapz(f(maskTot), p(maskTot));
Pin  = trapz(f(maskIn),  p(maskIn));
frac = Pin / Ptot;
end
