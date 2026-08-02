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
M  = simCfg.M; % 每控制周期内对象积分步数
h  = simCfg.h; % 对象积分步长

% 细网格: 共 N*M 步，将宿主结构的加速度输入的每一个变量乘以M倍，方便RK数值计算
xddFine = repelem(xddCtrl, M);   

% 预分配：控制器采样时刻的日志（长度 N）
log.t      = tCtrl;       % 时间 (s)
log.y      = zeros(N,1);  % DVA 相对位移 y
log.yd     = zeros(N,1);  % 相对速度
log.ydd    = zeros(N,1);  % 相对加速度（由对象方程算出，非数值微分）
log.ym     = zeros(N,1);  % 参考模型位移
log.ym_d   = zeros(N,1);  % 参考模型速度
log.ym_dd  = zeros(N,1);  % 参考模型加速度
log.e      = zeros(N,1);  % 位移跟踪误差 e = y - ym
log.ed     = zeros(N,1);  % 速度跟踪误差
log.edd    = zeros(N,1);  % 加速度跟踪误差
log.s      = zeros(N,1);  % 复合误差面 s = Bed + w1*Be + w2*a_f
log.a_f    = zeros(N,1);  % 经带通+低通后的带内加速度误差 a_f
log.u      = zeros(N,1);  % 主动控制力
log.dk_hat = zeros(N,1);  % 等效刚度增量估计 Δk_hat
log.dc_hat = zeros(N,1);  % 等效阻尼增量估计 Δc_hat
log.By     = zeros(N,1);  % 经带通+低通后的局部位移（自适应律用）
log.Byd    = zeros(N,1);  % 经带通+低通后的局部速度
log.xdd    = xddCtrl;     % 安装点/基座加速度激励（控制器采样网格）

% state：4x1 状态向量 [y, yd, ym, ym_d]
state = zeros(4,1);
% 初始化控制器状态 包括滤波器状态和自适应参数
stCtl = initControllerState(fltSys, ctl);
isPassive = strcmpi(ctl.mode, 'passive');
iFine = 1;

for k = 1:N
    % 测量时刻加速度：用上一拍保持的 u（首步为 0）
    if k == 1
        uHold = 0;
    else
        uHold = log.u(k-1);
    end
    xdd_k = xddCtrl(k);
    % 计算当前相对加速度和参考模型相对加速度，用于计算复合误差面
    % 真实实验中，ydd通过测量得到，ymdd通过测量得到
    [ydd, ymdd] = localPlantAccel(state, uHold, xdd_k, plant);

    % meas代表实该步可用的变量，其中一部分通过测量得到，参考模型的部分通过计算得到
    meas.y = state(1); meas.yd = state(2);
    meas.ym = state(3); meas.ym_d = state(4);
    meas.ydd = ydd; meas.ym_dd = ymdd;

    % rec：record的缩写，记录该步的控制器输出和状态
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
Fnl = plant.k1*y + plant.k3*y^3 + plant.k5*y^5;  % 刚度项提供的力，包括线性和非线性部分
ydd  = (-plant.c*yd  - Fnl - plant.m*xdd + u) / plant.m; % 根据现有状态计算相对加速度，在实际实验中通过测量得到
ymdd = (-plant.c_m*ymd - plant.k_m*ym - plant.m*xdd) / plant.m; % 参考模型相对加速度，通过参考模型方程计算得到
end
