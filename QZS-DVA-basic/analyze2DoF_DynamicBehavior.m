function analyze2DoF_DynamicBehavior(Omega_target, X_force_target, params)
% analyze2DoF_DynamicBehavior  二自由度 QZS-DVA 稳态非线性动力学诊断
%
%   analyze2DoF_DynamicBehavior(Omega_target, X_force_target)
%   analyze2DoF_DynamicBehavior(Omega_target, X_force_target, params)
%
%   通过 ode45 数值积分，在指定无量纲频率 Omega_target 与激励幅值
%   X_force_target 下，诊断系统的稳态时域波形、相轨迹及 FFT 频谱，
%   用于识别多倍频、次谐波或宽频混沌等非线性现象。
%
%   输入:
%       Omega_target   - 无量纲激励频率（相对 w1Residual）
%       X_force_target - 无量纲激励力幅值，覆盖 params.X_force
%       params         - （可选）由 getQzsDvaParams 返回的参数结构体
%
%   无量纲时间 tau = w1Residual * t_phys；
%   物理频率 f (Hz) = Omega * w1Residual / (2*pi)。

    if nargin < 3
        params = getQzsDvaParams('needDimensionless', true);
    end

    % 用目标激励幅值覆盖默认无量纲力
    params.X_force = X_force_target;

    w1Residual = params.w1Residual;

    %% ODE45 仿真配置
    numPeriods     = 500;  % 总仿真周期数，保证瞬态充分衰减
    extractPeriods = 100;   % 提取最后若干周期作为稳态段
    plotPeriods    = 30;   % 时域子图展示最后若干周期

    T_period = 2 * pi / Omega_target;           % 当前频率下一个激励周期（无量纲 tau）
    t_span   = [0, numPeriods * T_period];      % 总仿真时间区间
    Z_init   = [3; 0; 1.5; 0];                    % 全局静止启动

    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

    odeFuncs = FuncsQzsOde();
    [t, Z] = ode45(@(t, Z) odeFuncs.ode3rd(t, Z, Omega_target, params), ...
        t_span, Z_init, options);

    %% 稳态数据提取
    t_extract = (numPeriods - extractPeriods) * T_period;
    idxSteady = t >= t_extract;

    t_steady       = t(idxSteady);
    y_steady       = Z(idxSteady, 1);
    y_dot_steady   = Z(idxSteady, 2);
    x0_steady      = Z(idxSteady, 3);

    % 时域展示：最后 plotPeriods 个周期
    t_plot_start = (numPeriods - plotPeriods) * T_period;
    idxPlot = t >= t_plot_start;
    t_plot      = t(idxPlot);
    y_plot      = Z(idxPlot, 1);
    x0_plot     = Z(idxPlot, 3);

    % 物理时间轴（秒），便于横轴标注
    t_plot_phys = t_plot / w1Residual;

    %% FFT 频谱分析（针对局部位移 y_steady）
    N = numel(y_steady);
    dtau = mean(diff(t_steady));                % 无量纲时间步长（平均）
    fs_tau = 1 / dtau;                          % 以 tau 为自变量的采样率

    Y_fft = fft(y_steady);
    ampSpectrum = abs(Y_fft(1:floor(N/2))) * 2 / N;  % 单边幅值谱

    % 频率轴：先得到"每 tau 单位"的循环频率，再换算为物理 Hz
    % 关系: f_phys (Hz) = f_tau * w1Residual
    f_tau_axis = (0:floor(N/2)-1)' / N * fs_tau;
    f_axis     = f_tau_axis * w1Residual;

    % 基准激励物理频率（Hz）
    f_base = Omega_target * w1Residual / (2 * pi);

    %% 绘图：2x2 诊断面板
    figure('Name', '非线性动力学行为诊断', 'Color', 'w');

    % --- 子图 1：局部位移时域波形（最后 plotPeriods 个周期）---
    subplot(2, 2, 1);
    plot(t_plot_phys, y_plot, 'b-', 'LineWidth', 1.0);
    grid on;
    xlabel('物理时间 t (s)');
    ylabel('无量纲局部位移 y');
    title(sprintf('局部位移时域波形（最后 %d 个周期）', plotPeriods));

    % --- 子图 2：宿主位移时域波形 ---
    subplot(2, 2, 2);
    plot(t_plot_phys, x0_plot, 'r-', 'LineWidth', 1.0);
    grid on;
    xlabel('物理时间 t (s)');
    ylabel('无量纲宿主位移 \bar{x}_0');
    title(sprintf('宿主位移时域波形（最后 %d 个周期）', plotPeriods));

    % --- 子图 3：局部响应相轨迹 ---
    subplot(2, 2, 3);
    plot(y_steady, y_dot_steady, 'k-', 'LineWidth', 0.8);
    grid on;
    xlabel('无量纲局部位移 y');
    ylabel('无量纲局部位移速度 \dot{y}');
    title(sprintf('局部响应相轨迹（稳态 %d 个周期）', extractPeriods));
    axis equal;

    % --- 子图 4：局部响应 FFT 单边幅值谱 ---
    subplot(2, 2, 4);
    plot(f_axis, ampSpectrum, 'b-', 'LineWidth', 1.0);
    hold on;
    yMax = max(ampSpectrum);
    if yMax > 0
        plot([f_base, f_base], [0, yMax], 'r--', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('f_{base}=%.2f Hz', f_base));
    else
        xline(f_base, 'r--', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('f_{base}=%.2f Hz', f_base));
    end
    xlim([0, f_base*4]);
    hold off;
    grid on;
    xlabel('物理频率 f (Hz)');
    ylabel('单边幅值');
    title('局部位移 FFT 频谱（稳态段）');
    % legend('Location', 'best');

    % 总标题补充关键参数
    sgtitle(sprintf(['\\Omega=%.4f,  X_{force}=%.4f,  f_{base}=%.2f Hz'], ...
        Omega_target, X_force_target, f_base), 'FontWeight', 'normal');

    %% 命令行摘要
    A_y = (max(y_steady) - min(y_steady)) / 2;
    A_x0 = (max(x0_steady) - min(x0_steady)) / 2;
    fprintf('\n===== 非线性动力学诊断 =====\n');
    fprintf('Omega_target = %.4f,  X_force = %.4f\n', Omega_target, X_force_target);
    fprintf('f_base = %.4f Hz\n', f_base);
    fprintf('稳态局部位移幅值 A_y = %.6f\n', A_y);
    fprintf('稳态宿主位移幅值 A_x0 = %.6f\n', A_x0);
    fprintf('稳态周期数 = %d,  仿真总周期数 = %d\n', extractPeriods, numPeriods);
end
