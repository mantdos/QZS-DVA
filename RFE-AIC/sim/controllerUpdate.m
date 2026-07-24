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
    st     struct
    ctl    struct
    plant  struct
    fltSys struct
end

Ts = fltSys.Ts;
m  = plant.m;
km = plant.k_m;
cm = plant.c_m;

c1 = ctl.c1;
c2 = ctl.c2;
Ms = ctl.Ms;
K  = ctl.K;
wc = fltSys.omega_c;

% --- 7.1 误差量（加速度由方程给出，非数值微分）---
e   = meas.y  - meas.ym;
ed  = meas.yd - meas.ym_d;
edd = meas.ydd - meas.ym_dd;

% --- 7.2 带内提取（三路同构: BP → LP）---
[Be,  st.bp_e,  st.lp_e]  = localCascade(e,        st.bp_e,  st.lp_e,  fltSys);
[Bed, st.bp_ed, st.lp_ed] = localCascade(ed,       st.bp_ed, st.lp_ed, fltSys);
[a_f, st.bp_edd,st.lp_edd]= localCascade(edd,      st.bp_edd,st.lp_edd,fltSys);
[By,  st.bp_y,  st.lp_y]  = localCascade(meas.y,   st.bp_y,  st.lp_y,  fltSys);
[Byd, st.bp_yd, st.lp_yd] = localCascade(meas.yd,  st.bp_yd, st.lp_yd, fltSys);

% --- 7.3 复合误差面 ---
s = Bed + c1*Be + c2*a_f;

% --- 7.4 / 7.5 已知补偿 + 未知集总估计 ---
F_known   = -cm*Bed - km*Be;
F_unc_hat = st.dk_hat*By + st.dc_hat*Byd;

% --- 7.6 控制律 ---
u = (m/Ms) * (-K*s + c2*wc*a_f - c1*Bed) - F_known + F_unc_hat;

if isfinite(ctl.u_sat)
    u = max(-ctl.u_sat, min(ctl.u_sat, u));
end

% --- 7.7 自适应律（前向欧拉）---
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
% 带通 → 一阶低通
bp = fltSys.bp;
lp = fltSys.lp;
yBp = bp.Cd * xBp + bp.Dd * uIn;
xBp = bp.Ad * xBp + bp.Bd * uIn;
yOut = lp.Cd * xLp + lp.Dd * yBp;
xLp = lp.Ad * xLp + lp.Bd * yBp;
end
