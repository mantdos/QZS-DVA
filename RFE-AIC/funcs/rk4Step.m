function [stateNext, accel] = rk4Step(state, u, xdd, h, plant)
% rk4Step  单步四阶龙格-库塔积分（DVA + 参考模型）
%
%   [stateNext, accel] = rk4Step(state, u, xdd, h, plant)
%
%   状态: state = [y; yd; ym; ym_d]
%   对象:   m*ydd  + c*yd  + k1*y + k3*y^3 + k5*y^5 = -m*xdd + u
%   参考:   m*ymdd + c_m*ymd + k_m*ym                 = -m*xdd
%
%   accel 结构体返回当前步起点加速度: .ydd, .ymdd（用于误差计算）

arguments
    state (4,1) double
    u     (1,1) double
    xdd   (1,1) double
    h     (1,1) double {mustBePositive}
    plant struct
end

k1 = localF(state,            u, xdd, plant);
k2 = localF(state + 0.5*h*k1, u, xdd, plant);
k3 = localF(state + 0.5*h*k2, u, xdd, plant);
k4 = localF(state +     h*k3, u, xdd, plant);

stateNext = state + (h/6)*(k1 + 2*k2 + 2*k3 + k4);

% 步起点加速度（与控制更新时刻一致）
[ydd, ymdd] = localAccel(state, u, xdd, plant);
accel.ydd  = ydd;
accel.ymdd = ymdd;
end

%% ------------------------------------------------------------------------
function dX = localF(X, u, xdd, plant)
[ydd, ymdd] = localAccel(X, u, xdd, plant);
dX = [X(2); ydd; X(4); ymdd];
end

function [ydd, ymdd] = localAccel(X, u, xdd, plant)
y   = X(1);
yd  = X(2);
ym  = X(3);
ymd = X(4);

m  = plant.m;
c  = plant.c;
k1 = plant.k1;
k3 = plant.k3;
k5 = plant.k5;
km = plant.k_m;
cm = plant.c_m;

Fnl = k1*y + k3*y^3 + k5*y^5;
ydd  = (-c*yd  - Fnl      - m*xdd + u) / m;
ymdd = (-cm*ymd - km*ym   - m*xdd    ) / m;
end
