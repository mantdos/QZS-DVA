function log = simDvaOnly(plant, simCfg, fltSys, ctl, tCtrl, xddCtrl)
% simDvaOnly  DVA + 参考模型双时钟仿真
%
%   log = simDvaOnly(plant, simCfg, fltSys, ctl, tCtrl, xddCtrl)
%
%   双时钟:
%       外环: 控制器采样 fs_ctrl，u 零阶保持
%       内环: 每控制周期 M 步 RK4，步长 h = 1/(fs_ctrl*M)
%
%   ctl.mode = 'passive' 时强制 u=0（无控）

arguments
    plant   struct
    simCfg  struct
    fltSys  struct
    ctl     struct
    tCtrl   (:,1) double
    xddCtrl (:,1) double
end

N  = numel(tCtrl);
M  = simCfg.M;
h  = simCfg.h;

% 细网格: 共 N*M 步，覆盖 [t0, t0+N*Ts)
xddFine = repelem(xddCtrl, M);   % xdd 亦零阶保持（与 u 对齐，简单稳健）

% 预分配
log.t      = tCtrl;
log.y      = zeros(N,1);
log.yd     = zeros(N,1);
log.ydd    = zeros(N,1);
log.ym     = zeros(N,1);
log.ym_d   = zeros(N,1);
log.ym_dd  = zeros(N,1);
log.e      = zeros(N,1);
log.ed     = zeros(N,1);
log.edd    = zeros(N,1);
log.s      = zeros(N,1);
log.a_f    = zeros(N,1);
log.u      = zeros(N,1);
log.dk_hat = zeros(N,1);
log.dc_hat = zeros(N,1);
log.By     = zeros(N,1);
log.Byd    = zeros(N,1);
log.xdd    = xddCtrl;

state = zeros(4,1);
stCtl = initControllerState(fltSys, ctl);
isPassive = isfield(ctl, 'mode') && strcmpi(ctl.mode, 'passive');
iFine = 1;

for k = 1:N
    % 测量时刻加速度：用上一拍保持的 u（首步为 0）
    if k == 1
        uHold = 0;
    else
        uHold = log.u(k-1);
    end
    xdd_k = xddCtrl(k);
    [ydd, ymdd] = localPlantAccel(state, uHold, xdd_k, plant);

    meas.y = state(1); meas.yd = state(2);
    meas.ym = state(3); meas.ym_d = state(4);
    meas.ydd = ydd; meas.ym_dd = ymdd;

    if isPassive
        u_k = 0;
        rec.e   = meas.y  - meas.ym;
        rec.ed  = meas.yd - meas.ym_d;
        rec.edd = meas.ydd - meas.ym_dd;
        rec.s = 0; rec.a_f = 0; rec.By = 0; rec.Byd = 0;
    else
        [u_k, stCtl, rec] = controllerUpdate(meas, stCtl, ctl, plant, fltSys);
    end

    log.y(k)      = state(1);
    log.yd(k)     = state(2);
    log.ydd(k)    = ydd;
    log.ym(k)     = state(3);
    log.ym_d(k)   = state(4);
    log.ym_dd(k)  = ymdd;
    log.e(k)      = rec.e;
    log.ed(k)     = rec.ed;
    log.edd(k)    = rec.edd;
    log.s(k)      = rec.s;
    log.a_f(k)    = rec.a_f;
    log.u(k)      = u_k;
    log.dk_hat(k) = stCtl.dk_hat;
    log.dc_hat(k) = stCtl.dc_hat;
    log.By(k)     = rec.By;
    log.Byd(k)    = rec.Byd;

    for i = 1:M
        [state, ~] = rk4Step(state, u_k, xddFine(iFine), h, plant);
        iFine = iFine + 1;
    end
end
end

%% ------------------------------------------------------------------------
function [ydd, ymdd] = localPlantAccel(X, u, xdd, plant)
y = X(1); yd = X(2); ym = X(3); ymd = X(4);
Fnl = plant.k1*y + plant.k3*y^3 + plant.k5*y^5;
ydd  = (-plant.c*yd  - Fnl - plant.m*xdd + u) / plant.m;
ymdd = (-plant.c_m*ymd - plant.k_m*ym - plant.m*xdd) / plant.m;
end
