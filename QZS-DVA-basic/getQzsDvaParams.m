function params = getQzsDvaParams(options)
% getQzsDvaParams  返回 QZS-DVA 二自由度模型的统一物理参数
%
%   params = getQzsDvaParams()
%   params = getQzsDvaParams('isPrintSummary', true)
%
%   本函数为 QZS-DVA 文件夹内所有脚本/函数提供一致的模型常数，便于
%   幅频曲线复现及非线性对吸振性能影响的对比分析。
%
%   可选 Name-Value 参数:
%       isPrintSummary     - 是否在命令行打印吸振器固有频率摘要 (默认 false)
%       needDimensionless  - 是否追加 HBM 无量纲参数 (默认 false)
%       R_m                - 特征长度，用于位移无量纲化 (默认 0.001 m)
%       F_N                - 简谐激振力幅值 (默认 5 N)
%
%   输出 params 结构体字段:
%       m0, k0, c0, xi0  - 主结构质量、刚度、阻尼系数、阻尼比
%       m1, k1, c1, xi1  - 吸振器质量、线性正刚度、阻尼系数、阻尼比
%       kn1, kn3         - 吸振器非线性刚度一次项、三次项
%       w0, f0_Hz        - 主结构无阻尼固有圆频率及频率 (Hz)
%       w1Origin, f1Origin_Hz   - 吸振器原始线性固有频率
%       w1Residual, f1Residual_Hz - 残余刚度 (kn1) 对应的固有频率
%
%   当 needDimensionless 为 true 时，额外输出:
%       mu, lambda, alpha_n3, X_force, R_m, F_N
%
%   示例:
%       p = getQzsDvaParams();
%       p = getQzsDvaParams('isPrintSummary', true);
%       p = getQzsDvaParams('needDimensionless', true, 'R_m', 1e-3, 'F_N', 5);

arguments
    options.isPrintSummary (1, 1) logical = false
    options.needDimensionless (1, 1) logical = false
    options.R_m (1, 1) double {mustBePositive} = 1e-3
    options.F_N (1, 1) double {mustBePositive} = 5
end

% --- 主结构 ---
params.m0 = 10;        % 主结构质量 (kg)
params.k0 = 15000;     % 主结构刚度 (N/m)
params.xi0 = 0.003;    % 主结构阻尼比

% --- 吸振器 (QZS-DVA) ---
params.m1 = 2.5;       % 吸振器质量 (kg)
params.k1 = 3400;      % 吸振器线性正刚度 (N/m)
params.xi1 = 0.005;    % 吸振器阻尼比

% 吸振器非线性刚度: F_n = kn1 * x + kn3 * x^3
params.kn1 = 360.8;    % 非线性一次项刚度 (N/m)，准零刚度设计后的残余刚度
params.kn3 = 1.73e7;   % 非线性三次项刚度 (N/m^3)

% --- 由阻尼比导出的 viscous 阻尼系数 ---
params.c1 = 2 * params.xi1 * sqrt(params.k1 * params.m1);
params.c0 = 2 * params.xi0 * sqrt(params.k0 * params.m0);

% --- 便于幅频扫描的固有频率量 ---
params.w0 = sqrt(params.k0 / params.m0);
params.f0_Hz = params.w0 / (2 * pi);

params.w1Origin = sqrt(params.k1 / params.m1);
params.f1Origin_Hz = params.w1Origin / (2 * pi);

params.w1Residual = sqrt(params.kn1 / params.m1);
params.f1Residual_Hz = params.w1Residual / (2 * pi);

if options.needDimensionless
    % 以吸振器残余频率 w1Residual 为特征频率、R_m 为特征长度进行无量纲化
    params.R_m = options.R_m;
    params.F_N = options.F_N;

    params.mu = params.m1 / params.m0;
    params.lambda = params.w0 / params.w1Residual;
    params.alpha_n3 = params.kn3 * params.R_m^2 / params.kn1;
    params.alpha0 = params.k0 / params.kn1;
    params.X_force = params.F_N / (params.kn1 * params.R_m);
    % xi0、xi1 本身即为阻尼比，无需再换算
end

if options.isPrintSummary
    fprintf('吸振器原始固有频率: %.4f Hz\n', params.f1Origin_Hz);
    fprintf('吸振器残余刚度产生的固有频率: %.4f Hz\n', params.f1Residual_Hz);
    if options.needDimensionless
        fprintf('无量纲: mu=%.4f, lambda=%.4f, alpha_n3=%.4e, X_force=%.4f\n', ...
            params.mu, params.lambda, params.alpha_n3, params.X_force);
    end
end

end
