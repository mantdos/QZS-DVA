function fPeak = identifyResonanceFreq(plant, simCfg, fltSys, ctl, amp, opts)
% identifyResonanceFreq  简谐扫频识别 DVA 共振频率（≈宿主反共振）
%
%   fPeak = identifyResonanceFreq(plant, simCfg, fltSys, ctl, amp)
%   fPeak = identifyResonanceFreq(..., opts)
%
%   在 [fLo,fHi] 上扫频，以稳态 RMS(y) 峰值对应频率为共振频率。
%   控制律按 ctl 冻结/自适应/无控设定。

arguments
    plant  struct
    simCfg struct
    fltSys struct
    ctl    struct
    amp    (1,1) double
    opts.fLo     (1,1) double = NaN
    opts.fHi     (1,1) double = NaN
    opts.nFreq   (1,1) double = 15
    opts.T       (1,1) double = 25
    opts.tDiscard (1,1) double = 8
    opts.seed    (1,1) double = 1
end

f1 = plant.f1_Hz;
if isnan(opts.fLo), opts.fLo = max(0.8, 0.7*f1); end
if isnan(opts.fHi), opts.fHi = min(6.0, 2.2*f1); end

fGrid = linspace(opts.fLo, opts.fHi, opts.nFreq);
rmsY  = zeros(size(fGrid));

simLocal = simCfg;
simLocal.T = opts.T;
simLocal.t_discard = opts.tDiscard;

for i = 1:numel(fGrid)
    tCtrl = (0:1/simLocal.fs_ctrl:simLocal.T).';
    xdd = amp * sin(2*pi*fGrid(i)*tCtrl);
    log = simDvaOnly(plant, simLocal, fltSys, ctl, tCtrl, xdd);
    idx = log.t >= opts.tDiscard;
    rmsY(i) = rms(log.y(idx));
end

[~, iMax] = max(rmsY);
fPeak = fGrid(iMax);

% 抛物线细化（三点）
if iMax > 1 && iMax < numel(fGrid)
    f0 = fGrid(iMax-1); f1p = fGrid(iMax); f2 = fGrid(iMax+1);
    y0 = rmsY(iMax-1);  y1 = rmsY(iMax);   y2 = rmsY(iMax+1);
    denom = (y0 - 2*y1 + y2);
    if abs(denom) > eps
        fPeak = f1p + 0.5*(y0 - y2)/denom * (f2 - f0)/2;
    end
end
end
