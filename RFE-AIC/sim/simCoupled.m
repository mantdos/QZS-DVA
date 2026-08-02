function log = simCoupled(plant, host, simCfg, fltSys, ctl, tCtrl, Fctrl)
% simCoupled  宿主–DVA 耦合 + 参考模型双时钟仿真
%
%   log = simCoupled(plant, host, simCfg, fltSys, ctl, tCtrl, Fctrl)
%
%   双时钟（与 simDvaOnly 一致）:
%       外环: 控制器采样 fs_ctrl，u 零阶保持
%       内环: 每控制周期 M 步 RK4，步长 h = 1/(fs_ctrl*M)
%
%   输入激励为作用于宿主的力 F(t)，而非安装点加速度。
%   控制器接口不变：仍用 (y, yd, ydd, ym, ...) 调用 controllerUpdate。
%
%   ctl.mode = 'passive' 时强制 u=0。

arguments
    plant   struct
    host    struct
    simCfg  struct
    fltSys  struct
    ctl     struct
    tCtrl   (:,1) double
    Fctrl   (:,1) double
end

N  = numel(tCtrl);
M  = simCfg.M;
h  = simCfg.h;

Ffine = repelem(Fctrl, M);   % F 与 u 对齐零阶保持

log.t      = tCtrl;
log.x      = zeros(N,1);
log.xd     = zeros(N,1);
log.xdd    = zeros(N,1);
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
log.F      = Fctrl;
log.F_react = zeros(N,1);

state = zeros(6,1);
stCtl = initControllerState(fltSys, ctl);
isPassive = isfield(ctl, 'mode') && strcmpi(ctl.mode, 'passive');
iFine = 1;

for k = 1:N
    if k == 1
        uHold = 0;
    else
        uHold = log.u(k-1);
    end
    F_k = Fctrl(k);
    [xdd, ydd, ymdd, Freact] = localPlantAccel(state, uHold, F_k, plant, host);

    meas.y = state(3); meas.yd = state(4);
    meas.ym = state(5); meas.ym_d = state(6);
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

    log.x(k)       = state(1);
    log.xd(k)      = state(2);
    log.xdd(k)     = xdd;
    log.y(k)       = state(3);
    log.yd(k)      = state(4);
    log.ydd(k)     = ydd;
    log.ym(k)      = state(5);
    log.ym_d(k)    = state(6);
    log.ym_dd(k)   = ymdd;
    log.e(k)       = rec.e;
    log.ed(k)      = rec.ed;
    log.edd(k)     = rec.edd;
    log.s(k)       = rec.s;
    log.a_f(k)     = rec.a_f;
    log.u(k)       = u_k;
    log.dk_hat(k)  = stCtl.dk_hat;
    log.dc_hat(k)  = stCtl.dc_hat;
    log.By(k)      = rec.By;
    log.Byd(k)     = rec.Byd;
    log.F_react(k) = Freact;

    for i = 1:M
        [state, ~] = rk4StepCoupled(state, u_k, Ffine(iFine), h, plant, host);
        iFine = iFine + 1;
    end
end
end

%% ------------------------------------------------------------------------
function [xdd, ydd, ymdd, Freact] = localPlantAccel(X, u, F, plant, host)
x  = X(1); xd = X(2);
y  = X(3); yd = X(4);
ym = X(5); ymd = X(6);

f_q = plant.k1*y + plant.k3*y^3 + plant.k5*y^5;
Freact = plant.c*yd + f_q - u;
xdd  = (F - host.C*xd - host.K*x + Freact) / host.M;
ydd  = -xdd - (plant.c*yd + f_q - u) / plant.m;
ymdd = (-plant.c_m*ymd - plant.k_m*ym - plant.m*xdd) / plant.m;
end
