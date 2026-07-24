function [xdd, stats] = genExcitation(exc, t)
% genExcitation  生成窄带随机主导激励（基座加速度 xdd）
%
%   [xdd, stats] = genExcitation(exc, t)
%
%   输入:
%       exc - 激励参数结构体
%           .f_d        主导频率 Hz
%           .bw         窄带带宽 Hz（=0 且 rand_ratio=0 时退化为纯简谐）
%           .rand_ratio 随机分量功率占比 [0,1]
%           .amp        整体幅值标度（标定总 RMS）
%           .seed       随机种子
%       t   - 时间向量 (s)
%
%   输出:
%       xdd   - 基座加速度时程（DVA-only 方程中的 xdd）
%       stats - 实测统计: rms, psd_f, psd_p, A_h, A_r
%
%   说明:
%       第一版 DVA 单独仿真中，将激励直接作为基座加速度 xdd。

arguments
    exc struct
    t   (:,1) double
end

rng(exc.seed);

f_d = exc.f_d;
bw  = exc.bw;
r   = max(0, min(1, exc.rand_ratio)); % 随机分量功率占比 [0,1]；0=纯简谐，越大包络起伏越明显
amp = abs(exc.amp); % 总 RMS 幅值标度（m/s^2）
dt  = mean(diff(t));
fs  = 1 / dt;
n   = numel(t);

% --- 纯简谐退化（amp 标定 RMS，与混激励一致）---
if bw == 0 && r == 0
    xdd = amp * sqrt(2) * sin(2*pi*f_d*t);
    stats = localExcStats(xdd, fs, amp, 0);
    return;
end

% --- 单位 RMS 正弦 ---
sineRaw  = sin(2*pi*f_d*t);
sineUnit = sineRaw / max(rms(sineRaw), eps);

% --- 随机窄带：白噪声经二阶带通 ---
if r > 0
    if bw <= 0
        bwEff = max(0.02*f_d, 1e-3);   % 极窄带退化
    else
        bwEff = bw;
    end
    zeta = bwEff / (2*f_d);            % 带宽 ≈ 2*zeta*f_d
    zeta = max(zeta, 1e-4);
    wd   = 2*pi*f_d;
    % H(s) = 2*zeta*wd*s / (s^2 + 2*zeta*wd*s + wd^2)
    A = [0, 1; -wd^2, -2*zeta*wd];
    B = [0; 1];
    C = [0, 2*zeta*wd];
    D = 0;
    sys = ss(A, B, C, D);
    w = randn(n, 1);
    randRaw = lsim(sys, w, t);
    randUnit = randRaw / max(rms(randRaw), eps);
else
    randUnit = zeros(n, 1);
end

% 功率分配: P_rand/P_total = r
A_h = amp * sqrt(1 - r);
A_r = amp * sqrt(r);
xdd = A_h * sineUnit + A_r * randUnit;

stats = localExcStats(xdd, fs, A_h, A_r);
end

%% ------------------------------------------------------------------------
function stats = localExcStats(xdd, fs, A_h, A_r)
stats.rms = rms(xdd);
stats.A_h = A_h;
stats.A_r = A_r;
nSig = numel(xdd);
nWin = min(nSig, max(256, 2^floor(log2(max(nSig/4, 256)))));
nfft = max(256, 2^nextpow2(nWin));
if nWin >= nSig
    nWin = max(64, floor(nSig/2));
end
[psd_p, psd_f] = pwelch(xdd, hamming(nWin), round(nWin/2), nfft, fs);
stats.psd_f = psd_f;
stats.psd_p = psd_p;

fprintf('[genExcitation] RMS=%.4e, A_h=%.4e, A_r=%.4e\n', ...
    stats.rms, A_h, A_r);
end
