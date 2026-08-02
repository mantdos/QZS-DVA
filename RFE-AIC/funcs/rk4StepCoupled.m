function [stateNext, accel] = rk4StepCoupled(state, u, F, h, plant, host)
% rk4StepCoupled  耦合宿主+DVA+参考模型 单步 RK4
%
%   状态: [x; xd; y; yd; ym; ym_d]
%
%   显式耦合（无隐式代数环）:
%       f_q = k1*y + k3*y^3 + k5*y^5
%       xdd = ( F - C*xd - K*x + c*yd + f_q - u ) / M
%       ydd = -xdd - ( c*yd + f_q - u ) / m
%       参考模型由实测宿主加速度驱动:
%       ym_dd = (-c_m*ym_d - k_m*ym - m*xdd) / m

arguments
    state (6,1) double
    u     (1,1) double
    F     (1,1) double
    h     (1,1) double
    plant struct
    host  struct
end

k1 = localF(state,            u, F, plant, host);
k2 = localF(state + 0.5*h*k1, u, F, plant, host);
k3 = localF(state + 0.5*h*k2, u, F, plant, host);
k4 = localF(state +     h*k3, u, F, plant, host);

stateNext = state + (h/6)*(k1 + 2*k2 + 2*k3 + k4);

[xdd, ydd, ymdd, Freact] = localAccel(state, u, F, plant, host);
accel.xdd    = xdd;
accel.ydd    = ydd;
accel.ymdd   = ymdd;
accel.Freact = Freact;
end

%% ------------------------------------------------------------------------
function dX = localF(X, u, F, plant, host)
[xdd, ydd, ymdd] = localAccel(X, u, F, plant, host);
dX = [X(2); xdd; X(4); ydd; X(6); ymdd];
end

function [xdd, ydd, ymdd, Freact] = localAccel(X, u, F, plant, host)
x  = X(1); xd = X(2);
y  = X(3); yd = X(4);
ym = X(5); ymd = X(6);

m  = plant.m;
c  = plant.c;
k1 = plant.k1;
k3 = plant.k3;
k5 = plant.k5;
km = plant.k_m;
cm = plant.c_m;
M  = host.M;
C  = host.C;
K  = host.K;

f_q = k1*y + k3*y^3 + k5*y^5;
% F_react = c*yd + f_q - u  （与文档 2.5 节一致）
Freact = c*yd + f_q - u;
xdd  = (F - C*xd - K*x + Freact) / M;
ydd  = -xdd - (c*yd + f_q - u) / m;
ymdd = (-cm*ymd - km*ym - m*xdd) / m;
end
