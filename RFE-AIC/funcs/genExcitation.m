function [xdd, stats] = genExcitation(exc, t, host, flt)
% genExcitation  力经宿主传函 → 安装点加速度（单向传递）
%
%   [xdd, stats] = genExcitation(exc, t, host)
%   [xdd, stats] = genExcitation(exc, t, host, flt)
%
%   物理背景:
%       实际激励是作用在宿主结构上的力，DVA 感受到的是安装点加速度。
%       安装点加速度的窄带特性来自宿主在 f_n 附近的共振放大，而非人为设定带宽。
%
%       本阶段采用单向传递（宿主 → DVA），忽略 DVA 反作用力对宿主的影响。
%       对应质量比较小的情形，用于控制器验证。
%       注意：单向传递下 xdd 不含反共振凹口，因此本阶段不能评估减振量，
%       只能评估跟踪性能（e、edd、参数收敛）。减振量需等耦合仿真。
%
%   生成流程:
%       F(t) = A_h*sin(2*pi*f_d*t) + F_r(t)   % 力端功率比 = rand_ratio
%         ↓  H(s) = s^2 / (M*s^2 + C*s + K)   % 力 → 绝对加速度
%       xdd_raw(t)
%         ↓  归一化（必须在宿主滤波之后）
%       xdd = xdd_raw / rms(xdd_raw) * exc.amp
%
%   输入:
%       exc - 激励参数
%           .f_d        谐波频率 Hz
%           .amp        安装点加速度 RMS 目标值 (m/s^2)
%           .rand_ratio 力端随机功率占比 [0,1]
%           .f_lo       随机分量带限下界 Hz（防 DC 推离 QZS 工作点）
%           .f_hi       随机分量带限上界 Hz
%           .seed       随机种子
%       t    - 时间向量 (s)
%       host - 宿主参数（M, C, K；可由 buildHost 得到）
%       flt  - （可选）控制器滤波参数，用于带内诊断 r_inband
%
%   输出:
%       xdd   - 安装点加速度时程 (m/s^2)，喂给 simDvaOnly
%       stats - rms / psd_f / psd_p / F_scale / r_inband
%
%   为何必须带限、不能用裸 randn:
%       1) randn 功率平摊在 [0, fs/2]；fs=2000 时落在 f_d 附近占比极小，
%          使 rand_ratio 失去物理意义，且结果随采样率变化。
%       2) 高通 f_lo 必须有：QZS 在平衡点附近刚度接近零，DC/准静态分量
%          会把工作点推离低刚度区，仿真反映的是工作点偏移而非 k3/k5 非线性。
%
%   滤波器阶次说明（避免与控制器混淆）:
%       力端随机用 4 阶 Butterworth 追求陡峭滚降，无相对阶约束。
%       控制器内 designBandpass 必须相对阶为 1（推导要求），二者职责不同。

arguments
    exc  struct
    t    (:,1) double
    host struct
    flt  struct = struct()
end

rng(exc.seed);

f_d = exc.f_d;
r   = max(0, min(1, exc.rand_ratio));
amp = abs(exc.amp);
dt  = mean(diff(t));
fs  = 1 / dt;
n   = numel(t);

% --- 力端：单位 RMS 谐波 + 带限随机，按 rand_ratio 分配功率 ---
sineUnit = sin(2*pi*f_d*t);
sineUnit = sineUnit / rms(sineUnit);

if r == 0
    F = sineUnit;
    Fr = zeros(n, 1);
else
    % --- 带限宽带随机：4 阶 Butterworth [f_lo, f_hi] ---
    % SOS 避免 b/a 在 fs>>f_hi 时病态；filtfilt 零相位，避免瞬态
    f_lo = exc.f_lo;
    f_hi = exc.f_hi;
    nyq  = fs / 2;
    [z, p, k] = butter(4, [f_lo, f_hi] / nyq, 'bandpass');
    [sos, g] = zp2sos(z, p, k);
    randRaw = filtfilt(sos, g, randn(n, 1));
    randUnit = randRaw / rms(randRaw);
    % 功率分配: P_rand/P_total = r（力端）
    A_h = sqrt(1 - r);   % 谐波分量单位-RMS 权重
    A_r = sqrt(r);       % 随机分量单位-RMS 权重
    F  = A_h * sineUnit + A_r * randUnit;
    Fr = A_r * randUnit; % 随机分量时程（单独过宿主，用于 r_inband）
end

% --- 宿主传函：力 → 绝对加速度（状态空间 + lsim）---
% M*xdd + C*xd + K*x = F  →  xdd = (F - C*xd - K*x)/M
% H(s) = s^2 / (M*s^2 + C*s + K)，D≠0，高频直通，物理正确
Ah = [0, 1; -host.K/host.M, -host.C/host.M];
Bh = [0; 1/host.M];
Ch = [-host.K/host.M, -host.C/host.M];
Dh = 1/host.M;
sysHost = ss(Ah, Bh, Ch, Dh);

% 激励过宿主结构，获得宿主结构加速度激励
xddRaw = lsim(sysHost, F, t);
xddRaw = xddRaw(:);

% 随机分量激励过宿主结构，获得宿主结构加速度激励
xddRandRaw = lsim(sysHost, Fr, t);
xddRandRaw = xddRandRaw(:);

% --- 输出端归一化：exc.amp = 安装点加速度 RMS ---
rmsRaw = rms(xddRaw);
if rmsRaw < eps
    error('genExcitation: xdd_raw RMS 过小，检查宿主参数或激励频率');
end
F_scale = amp / rmsRaw;
xdd     = xddRaw * F_scale; % 宿主结构加速度激励归一化
xddRand = xddRandRaw * F_scale; % 宿主结构加速度激励归一化（随机分量）

stats = localExcStats(xdd, xddRand, fs, f_d, F_scale, flt);
end

%% ------------------------------------------------------------------------
function stats = localExcStats(xdd, xddRand, fs, f_d, F_scale, flt)
stats.rms     = rms(xdd);
stats.F_scale = F_scale;

nSig = numel(xdd);
nWin = min(nSig, max(256, 2^floor(log2(max(nSig/4, 256)))));
nfft = max(256, 2^nextpow2(nWin));
if nWin >= nSig
    nWin = max(64, floor(nSig/2));
end
win = hamming(nWin);
nover = round(nWin/2);
[psd_p, psd_f] = pwelch(xdd, win, nover, nfft, fs);
stats.psd_f = psd_f; % 对应的频率
stats.psd_p = psd_p; % PSD，单位为信号单位^2/Hz

% 带内扰动诊断：论文工况应报告 r_inband，而非力端总功率占比 rand_ratio
% 宿主放大后谐波更突出，r_inband 通常小于 rand_ratio
stats.r_inband = NaN;
% 带通滤波器半带宽：Delta = zeta_b*omega_d/(2*pi)
dF = flt.zeta_b * flt.omega_d / (2*pi);
f1 = f_d - dF;
f2 = f_d + dF;
[Pr, fr] = pwelch(xddRand, win, nover, nfft, fs);
pRand = localBandPower(fr, Pr, f1, f2); % [f1,f2] 上对随机分量 PSD 积分
pTot = localBandPower(psd_f, psd_p, f1, f2); % [f1,f2] 上对总信号 PSD 积分
if pTot > 0
    stats.r_inband = pRand / pTot;
end

fprintf('[genExcitation] rms(xdd)=%.4e, F_scale=%.4e, r_inband=%.4f\n', ...
    stats.rms, stats.F_scale, stats.r_inband);
end

function P = localBandPower(f, pxx, f1, f2)
mask = f >= f1 & f <= f2;
if ~any(mask)
    P = 0;
else
    P = trapz(f(mask), pxx(mask));
end
end
