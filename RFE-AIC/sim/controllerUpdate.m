function [u, st, rec] = controllerUpdate(meas, st, ctl, plant, fltSys)
% controllerUpdate  单步离散控制器更新（滤波 + 控制律 + 自适应）
%
%   [u, st, rec] = controllerUpdate(meas, st, ctl, plant, fltSys)
%
%   meas 字段: y, yd, ydd, ym, ym_d, ym_dd
%   st   维护: 五路滤波器状态 + dk_hat + dc_hat
%
%   第一版: 裸自适应律，无投影、无 σ-modification、无死区。

arguments
    meas   struct
    st     struct   % 所有五路通道的滤波器状态变量，以及两个自适应参数
    ctl    struct
    plant  struct
    fltSys struct
end

Ts = fltSys.Ts;
m  = plant.m;
km = plant.k_m;
cm = plant.c_m;

w1 = ctl.w1;
w2 = ctl.w2;
Ms = ctl.Ms;
K  = ctl.K;
wc = fltSys.omega_c;

% --- 7.1 误差量（加速度由方程给出，非数值微分）---
e   = meas.y  - meas.ym;
ed  = meas.yd - meas.ym_d;
edd = meas.ydd - meas.ym_dd;

% --- 7.2 、e、ed、y、yd：只过带通
[Be,  st.bp_e]  = localBp(e,        st.bp_e,  fltSys);
[Bed, st.bp_ed] = localBp(ed,       st.bp_ed, fltSys);
[By,  st.bp_y]  = localBp(meas.y,   st.bp_y,  fltSys);
[Byd, st.bp_yd] = localBp(meas.yd,  st.bp_yd, fltSys);
% edd：带通 → 一阶低通
[a_f_raw, st.lp_edd] = localLp(edd,     st.lp_edd, fltSys);   % 控制律用
[a_f,     st.bp_edd] = localBp(a_f_raw, st.bp_edd, fltSys);   % 自适应律用


% --- 7.3 / 7.4 / 7.5 误差面 + 已知补偿 + 未知集总估计 ---
% 全部使用未滤波量，因为控制力u直接作用于动力学系统，不需要经过滤波器
s_raw     = ed + w1*e + w2*a_f_raw;
F_known   = -cm*ed - km*e;
F_unc_hat = st.dk_hat*meas.y + st.dc_hat*meas.yd;

% --- 7.6 控制律 全部使用未滤波量 ---
u = (m/Ms) * (-K*s_raw + w2*wc*a_f_raw - w1*ed) - F_known + F_unc_hat;

% u = 0;

% --- 7.7 自适应律（前向欧拉）---
s = Bed + w1*Be + w2*a_f;
if ctl.adapt_on
    st.dk_hat = st.dk_hat - ctl.Gamma_k * s * By  * Ts;
    st.dc_hat = st.dc_hat - ctl.Gamma_c * s * Byd * Ts;
else
    st.dk_hat = ctl.dk_frozen;
    st.dc_hat = ctl.dc_frozen;
end

% --- 7.8 记录 ---
rec.e = e; rec.ed = ed; rec.edd = edd;
rec.s = s; rec.a_f = a_f;
rec.Be = Be; rec.Bed = Bed;
rec.By = By; rec.Byd = Byd;
rec.u = u;
rec.dk_hat = st.dk_hat;
rec.dc_hat = st.dc_hat;
rec.F_known = F_known;
rec.F_unc_hat = F_unc_hat;
end

%% ------------------------------------------------------------------------
function [yOut, xBp, xLp] = localCascade(uIn, xBp, xLp, fltSys)
%  离散状态空间标准形式：x[k+1]=Ad?x[k]+Bd?u[k]，y[k]=Cd?x[k]+Dd?u[k]
% 带通 → 一阶低通
bp = fltSys.bp;
lp = fltSys.lp;
yBp = bp.Cd * xBp + bp.Dd * uIn; % 输出方程：y[k] = Cd·x[k] + Dd·u[k]，根据当前时刻的输入和状态，计算当前时刻的输出
xBp = bp.Ad * xBp + bp.Bd * uIn; % 状态方程：x[k+1] = Ad·x[k] + Bd·u[k]，根据当前时刻的输入和状态，提前计算下一时刻的状态，并回写到xBp中
yOut = lp.Cd * xLp + lp.Dd * yBp;
xLp = lp.Ad * xLp + lp.Bd * yBp;
end

function [yOut, xLp] = localLp(uIn, xLp, fltSys)
    %  一阶低通
    lp = fltSys.lp;
    yOut = lp.Cd * xLp + lp.Dd * uIn;
    xLp = lp.Ad * xLp + lp.Bd * uIn;
    end


function [yOut, xBp] = localBp(uIn, xBp, fltSys)
bp = fltSys.bp;
yOut = bp.Cd * xBp + bp.Dd * uIn;
xBp = bp.Ad * xBp + bp.Bd * uIn;
end