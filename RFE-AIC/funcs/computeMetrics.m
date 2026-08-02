function metrics = computeMetrics(log, plant, simCfg)
% computeMetrics  跟踪性能指标（单向传递验证阶段）
%
%   metrics = computeMetrics(log, plant, simCfg)
%
%   保留:
%       e_rms, ed_rms, edd_rms, y_rms, u_rms
%       u_over_Fpas  — ||u||/||f_qzs||，验证主动力只做微调
%       dk_mean      — 稳态段 dk_hat 均值
%
%   单向传递下不计算减振量 / 反共振等耦合指标。

arguments
    log    struct
    plant  struct
    simCfg struct
end

t   = log.t(:);
idx = t >= simCfg.t_discard;
y   = log.y(idx);
e   = log.e(idx);
ed  = log.ed(idx);
edd = log.edd(idx);
u   = log.u(idx);

metrics.y_rms   = rms(y);
metrics.e_rms   = rms(e);
metrics.ed_rms  = rms(ed);
metrics.edd_rms = rms(edd);
metrics.u_rms   = rms(u);

Fpas = plant.k1*y + plant.k3*y.^3 + plant.k5*y.^5;
metrics.u_over_Fpas = metrics.u_rms / max(rms(Fpas), eps);
metrics.dk_mean = mean(log.dk_hat(idx));
end
