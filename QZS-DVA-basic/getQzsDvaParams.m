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
%       isPrintSummary - 是否在命令行打印吸振器固有频率摘要 (默认 false)
%
%   输出 params 结构体字段:
%       m0, k0, c0, xi0  - 主结构质量、刚度、阻尼系数、阻尼比
%       m1, k1, c1, xi1  - 吸振器质量、线性正刚度、阻尼系数、阻尼比
%       kn1, kn3         - 吸振器非线性刚度一次项、三次项
%       w0, f0_Hz        - 主结构无阻尼固有圆频率及频率 (Hz)
%       w1Origin, f1Origin_Hz   - 吸振器原始线性固有频率
%       w1Residual, f1Residual_Hz - 残余刚度 (kn1) 对应的固有频率
%
%   示例:
%       p = getQzsDvaParams();
%       p = getQzsDvaParams('isPrintSummary', true);

arguments
    options.isPrintSummary (1, 1) logical = false
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

if options.isPrintSummary
    fprintf('吸振器原始固有频率: %.4f Hz\n', params.f1Origin_Hz);
    fprintf('吸振器残余刚度产生的固有频率: %.4f Hz\n', params.f1Residual_Hz);
end

end
