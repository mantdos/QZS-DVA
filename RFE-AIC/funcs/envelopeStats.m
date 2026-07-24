function stats = envelopeStats(y, t, tDiscard)
% envelopeStats  信号包络统计诊断
%
%   stats = envelopeStats(y, t, tDiscard)
%
%   输出字段:
%       A_mean   - 包络均值
%       A_std    - 包络标准差
%       A_cv     - 变异系数 σ_A / ?
%       tau_env  - 包络自相关衰减到 1/e 的时间 (s)
%       f_env_bw - 包络功率谱等效带宽 (Hz)
%
%   判据提示:
%       A_cv < 10%  → 主导正弦占优，等效刚度近似定值
%       A_cv > 30%  → 随机显著，存在固有瞬时失配

arguments
    y         (:,1) double
    t         (:,1) double
    tDiscard  (1,1) double {mustBeNonnegative} = 0
end

idx = t >= tDiscard;
y = y(idx);
t = t(idx);
dt = mean(diff(t));
fs = 1 / dt;

% Hilbert 包络
A = abs(hilbert(y));
A = A(:);

stats.A_mean = mean(A);
stats.A_std  = std(A);
stats.A_cv   = stats.A_std / max(stats.A_mean, eps);

% 包络自相关 → 1/e 时间
A0 = A - mean(A);
maxLag = min(numel(A0)-1, round(fs * min(30, 0.25*numel(A0)/fs)));
[ac, lags] = xcorr(A0, maxLag, 'coeff');
acPos = ac(lags >= 0);
tLag  = lags(lags >= 0) / fs;
thr   = exp(-1);
iCross = find(acPos < thr, 1, 'first');
if isempty(iCross) || iCross < 2
    stats.tau_env = NaN;
else
    % 线性插值过 1/e
    a1 = acPos(iCross-1); a2 = acPos(iCross);
    t1 = tLag(iCross-1);  t2 = tLag(iCross);
    stats.tau_env = t1 + (thr-a1)/(a2-a1+eps)*(t2-t1);
end

% 包络功率谱等效带宽: (∫G df)^2 / ∫G^2 df
nSig = numel(A);
nWin = min(nSig, max(256, 2^floor(log2(max(nSig/4, 256)))));
if nWin >= nSig
    nWin = max(64, floor(nSig/2));
end
nfft = max(256, 2^nextpow2(nWin));
[Gp, fp] = pwelch(A0, hamming(nWin), round(nWin/2), nfft, fs);
df = mean(diff(fp));
Gsum  = sum(Gp) * df;
G2sum = sum(Gp.^2) * df;
stats.f_env_bw = (Gsum^2) / max(G2sum, eps);

fprintf(['[envelopeStats] A_mean=%.4e, A_cv=%.2f%%, tau_env=%.3f s, ' ...
    'f_env_bw=%.4f Hz\n'], stats.A_mean, 100*stats.A_cv, ...
    stats.tau_env, stats.f_env_bw);

if stats.A_cv < 0.10
    fprintf('  → 主导正弦占优，等效刚度近似定值\n');
elseif stats.A_cv > 0.30
    fprintf('  → 随机显著，存在固有瞬时失配\n');
else
    fprintf('  → 混合激励区间\n');
end
end
