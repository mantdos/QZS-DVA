function fltSys = designBandpass(flt, fsCtrl)
% designBandpass  设计二阶带通 + 一阶低通（Tustin 预畸变）
%
%   fltSys = designBandpass(flt, fsCtrl)
%
%   连续带通（相对阶 = 1）:
%       H_bp(s) = 2*zeta_b*omega_d*s / (s^2 + 2*zeta_b*omega_d*s + omega_d^2)
%   连续低通:
%       H_lp(s) = omega_c / (s + omega_c)
%
%   输出 fltSys 字段:
%       .bp  - 带通离散 SS: Ad,Bd,Cd,Dd, nx
%       .lp  - 低通离散 SS: Ad,Bd,Cd,Dd, nx
%       .phaseLag_deg - 级联在 omega_d 处的相位滞后 (°)
%
%   三路同构: e / ? / ? 各自独立实例，均为 带通 → 一阶低通。

arguments
    flt    struct
    fsCtrl (1,1) double
end

wd = flt.omega_d;
zb = flt.zeta_b;
wc = flt.omega_c;
Ts = 1 / fsCtrl;

% --- 连续带通可控标准型 ---
Abp = [0, 1; -wd^2, -2*zb*wd];
Bbp = [0; 1];
Cbp = [0, 2*zb*wd];
Dbp = 0;

% --- 连续一阶低通 ---
Alp = -wc;
Blp = wc;   % ? = -wc*x + wc*u, y = x  → H = wc/(s+wc)
Clp = 1;
Dlp = 0;

bp = localC2dTustin(Abp, Bbp, Cbp, Dbp, Ts, wd);
lp = localC2dTustin(Alp, Blp, Clp, Dlp, Ts, wc);

fltSys.bp = bp;
fltSys.lp = lp;
fltSys.Ts = Ts;
fltSys.omega_d = wd;
fltSys.omega_c = wc;
fltSys.zeta_b = zb;

% --- 级联在 omega_d 处相位滞后 ---
z = exp(1j*wd*Ts);
Hbp = bp.Cd * ((z*eye(bp.nx) - bp.Ad) \ bp.Bd) + bp.Dd;
Hlp = lp.Cd * ((z*eye(lp.nx) - lp.Ad) \ lp.Bd) + lp.Dd;
Htot = Hbp * Hlp;
phaseLag = -rad2deg(angle(Htot));
fltSys.phaseLag_deg = phaseLag;
fltSys.mag_at_wd = abs(Htot);

fprintf('[designBandpass] omega_d=%.3f rad/s (%.3f Hz), zeta_b=%.3f, omega_c=%.3f\n', ...
    wd, wd/(2*pi), zb, wc);
fprintf('[designBandpass] 级联在 omega_d 处: |H|=%.4f, 相位滞后=%.2f°', ...
    abs(Htot), phaseLag);
if phaseLag < 20
    fprintf('  (<20° OK)\n');
else
    fprintf('  (!! 超过 20°，建议增大 omega_c 或减小 zeta_b)\n');
end
end

%% ------------------------------------------------------------------------
function sysd = localC2dTustin(A, B, C, D, Ts, wPre)
sysc = ss(A, B, C, D);
op = c2dOptions('Method', 'tustin', 'PrewarpFrequency', wPre);
sysdd = c2d(sysc, Ts, op);
sysd.Ad = sysdd.A;
sysd.Bd = sysdd.B;
sysd.Cd = sysdd.C;
sysd.Dd = sysdd.D;
sysd.nx = size(sysdd.A, 1);
end
