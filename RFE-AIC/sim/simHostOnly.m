function log = simHostOnly(host, simCfg, tCtrl, Fctrl)
% simHostOnly  无 DVA 宿主基线响应
%
%   log = simHostOnly(host, simCfg, tCtrl, Fctrl)
%
%   方程: M*xdd + C*xd + K*x = F(t)
%   双时钟细分 RK4（与耦合仿真相同步长），用于计算减振量分母。

arguments
    host   struct
    simCfg struct
    tCtrl  (:,1) double
    Fctrl  (:,1) double
end

N = numel(tCtrl);
M = simCfg.M;
h = simCfg.h;
Ffine = repelem(Fctrl, M);

log.t   = tCtrl;
log.x   = zeros(N,1);
log.xd  = zeros(N,1);
log.xdd = zeros(N,1);
log.F   = Fctrl;

state = zeros(2,1);
iFine = 1;

for k = 1:N
    F_k = Fctrl(k);
    xdd = (F_k - host.C*state(2) - host.K*state(1)) / host.M;
    log.x(k)   = state(1);
    log.xd(k)  = state(2);
    log.xdd(k) = xdd;

    for i = 1:M
        state = localRk4(state, Ffine(iFine), h, host);
        iFine = iFine + 1;
    end
end
end

%% ------------------------------------------------------------------------
function stateNext = localRk4(state, F, h, host)
k1 = localF(state,            F, host);
k2 = localF(state + 0.5*h*k1, F, host);
k3 = localF(state + 0.5*h*k2, F, host);
k4 = localF(state +     h*k3, F, host);
stateNext = state + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function dX = localF(X, F, host)
xdd = (F - host.C*X(2) - host.K*X(1)) / host.M;
dX = [X(2); xdd];
end
