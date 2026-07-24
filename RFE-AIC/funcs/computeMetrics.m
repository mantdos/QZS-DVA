function metrics = computeMetrics(log, plant, simCfg)
% computeMetrics  统一性能指标（所有实验共用）
%
%   metrics = computeMetrics(log, plant, simCfg)
%
%   指标:
%       f_ar_Hz      - 反共振/等效固有频率估计 (Hz)（由响应谱峰）
%       y_rms        - 吸振器相对位移 RMS
%       e_rms        - ‖e‖_rms
%       edd_rms      - ‖?‖_rms
%       u_rms        - 控制力 RMS
%       Fpas_rms     - 被动恢复力 RMS（k1*y+k3*y^3+k5*y^5）
%       u_over_Fpas  - u_rms / Fpas_rms
%       dk_mean      - 稳态段 dk_hat 均值
%       dc_mean      - 稳态段 dc_hat 均值

arguments
    log     struct
    plant   struct
    simCfg  struct
end

t  = log.t(:);
idx = t >= simCfg.t_discard;
t2  = t(idx);
y   = log.y(idx);
e   = log.e(idx);
edd = log.edd(idx);
u   = log.u(idx);
fs  = 1 / mean(diff(t2));

metrics.y_rms   = rms(y);
metrics.e_rms   = rms(e);
metrics.edd_rms = rms(edd);
metrics.u_rms   = rms(u);

Fpas = plant.k1*y + plant.k3*y.^3 + plant.k5*y.^5;
metrics.Fpas_rms = rms(Fpas);
metrics.u_over_Fpas = metrics.u_rms / max(metrics.Fpas_rms, eps);

if isfield(log, 'dk_hat')
    metrics.dk_mean = mean(log.dk_hat(idx));
    metrics.dc_mean = mean(log.dc_hat(idx));
else
    metrics.dk_mean = NaN;
    metrics.dc_mean = NaN;
end

% 响应谱峰 → 等效固有/反共振频率估计
nSig = numel(y);
nWin = min(nSig, max(256, 2^floor(log2(max(nSig/4, 256)))));
if nWin >= nSig
    nWin = max(64, floor(nSig/2));
end
nfft = max(256, 2^nextpow2(nWin));
[Py, fy] = pwelch(y - mean(y), hamming(nWin), round(nWin/2), nfft, fs);
% 限制在合理频带内搜索
fLo = max(0.5, plant.f1_Hz * 0.5);
fHi = plant.f1_Hz * 2.5;
mask = fy >= fLo & fy <= fHi;
if any(mask)
    [~, iMax] = max(Py(mask));
    fyM = fy(mask);
    metrics.f_ar_Hz = fyM(iMax);
else
    [~, iMax] = max(Py);
    metrics.f_ar_Hz = fy(iMax);
end

% 由包络幅值给出的描述函数固有频率（辅助）
A = mean(abs(hilbert(y)));
dkA = (3/4)*plant.k3*A^2 + (5/8)*plant.k5*A^4;
metrics.f_df_Hz = sqrt((plant.k1 + dkA) / plant.m) / (2*pi);
metrics.A_env = A;
metrics.sigma_y = std(y);
end
