function Funcs = FuncsQzsOde()
% FuncsQzsOde  二自由度 QZS-DVA 无量纲状态空间 ODE 工具包
%
%   Funcs = FuncsQzsOde()
%
%   返回结构体 Funcs，包含不同非线性刚度阶次下的动力学方程句柄。
%   自变量为无量纲时间 tau（tau = w1Residual * t_phys），状态向量
%   Z = [y; y_dot; x0_bar; x0_bar_dot]。

    % 纯立方刚度（三阶）非线性 QZS-DVA 动力学方程
    Funcs.ode3rd = @ode3rd;
    % 包含五阶刚度的非线性 QZS-DVA 动力学方程（预留）
    Funcs.ode5th = @ode5th;
end

function dZdt = ode3rd(t, Z, Omega, params)
% ode3rd  三阶（立方）非线性 QZS-DVA 二自由度状态方程
%
%   dZdt = ode3rd(t, Z, Omega, params)
%
%   输入:
%       t      - 无量纲时间 tau
%       Z      - 状态向量 [y; y_dot; x0_bar; x0_bar_dot]
%       Omega  - 无量纲激励频率（相对 w1Residual）
%       params - 含 mu, xi1, xi0, lambda, alpha0, alpha_n3, X_force 的结构体
%
%   方程（与一阶 HBM 推导一致）:
%       dZ1 = Z2
%       dZ3 = Z4
%       dZ4 = 2*mu*xi1*Z2 + mu*Z1 + mu*alpha_n3*Z1^3 ...
%             - 2*xi0*lambda*Z4 - mu*alpha0*Z3 + X_force*sin(Omega*t)
%       dZ2 = -dZ4 - 2*xi1*Z2 - Z1 - alpha_n3*Z1^3

    mu       = params.mu;
    xi1      = params.xi1;
    xi0      = params.xi0;
    lambda   = params.lambda;
    alpha0   = params.alpha0;
    alpha_n3 = params.alpha_n3;
    X_force  = params.X_force;

    Z1 = Z(1);
    Z2 = Z(2);
    Z3 = Z(3);
    Z4 = Z(4);

    dZ4 = 2 * mu * xi1 * Z2 + mu * Z1 + mu * alpha_n3 * Z1^3 ...
        - 2 * xi0 * lambda * Z4 - mu * alpha0 * Z3 + X_force * sin(Omega * t);
    dZ2 = -dZ4 - 2 * xi1 * Z2 - Z1 - alpha_n3 * Z1^3;

    dZdt = [Z2; dZ2; Z4; dZ4];
end

function dZdt = ode5th(t, Z, Omega, params)
% ode5th  含五阶刚度的 QZS-DVA 方程（占位，待后续推导补充）

    warning('5阶刚度方程待后续推导补充');
    dZdt = zeros(4, 1);
end
