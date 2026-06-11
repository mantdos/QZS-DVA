function Funcs = FuncsDVAUtil()
    %% 收拢DVA二自由度系统的结构参数信息，并生成相应的函数，防止每个mat文件都要重复输入
    Funcs.getDVAAndMainSystemParams = @getDVAAndMainSystemParams;

    %% 传递函数与状态转移矩阵
    % 获取二自由度DVA的f0→a_0,f0→a_1的传递函数(加速度)
    Funcs.calcuTransferAccAmplitude = @calcuTransferAccAmplitude;
    % 获取二自由度DVA的f0→v_0,f0→v_1的传递函数(速度)
    Funcs.calcuTransferVelocityAmplitude = @calcuTransferVelocityAmplitude;
    % 获取二自由度DVA状态方程的状态矩阵和输入矩阵
    Funcs.getDVAWithMainSystemStateMatrix = @getDVAWithMainSystemStateMatrix;

    %% 谱能量相关
    % 获取二自由度DVA的x0'→x1'的谱能量比指标(解析解)
    Funcs.calcuSpectralEnergyRatioUnderRandomExcite = @calcuSpectralEnergyRatioUnderRandomExcite;
    % 获取二自由度DVA的x0'→x1'的谱能量比指标关于k1的极大值点(解析解)
    Funcs.calcuMaximumK1SpectralEnergyRatioUnderRandomExcite = @calcuMaximumK1SpectralEnergyRatioUnderRandomExcite;
    % 获取二自由度DVA的f0到x0',x1'的谱能量指标(H2指标)(解析解)
    Funcs.calcuSpectralEnergyUnderRandomExcite = @calcuSpectralEnergyUnderRandomExcite;

     % 获取二自由度DVA的x0→x1的谱能量比指标(解析解)
     Funcs.calcuSpectralEnergyRatioUnderRandomExciteWithDis = @calcuSpectralEnergyRatioUnderRandomExciteWithDis;
     % 获取二自由度DVA的x0→x1的谱能量比指标关于k1的极大值点(解析解)
     Funcs.calcuMaximumK1SpectralEnergyRatioUnderRandomExciteWithDis = @calcuMaximumK1SpectralEnergyRatioUnderRandomExciteWithDis;
     % 获取二自由度DVA的f0到x0,x1的谱能量指标(H2指标)(解析解)
     Funcs.calcuSpectralEnergyUnderRandomExciteWithDis = @calcuSpectralEnergyUnderRandomExciteWithDis;

    % 根据两个时域信号计算其各自RMS值及谱能量比
    Funcs.calcRMSAndSpectrumEnergyRatio = @calcRMSAndSpectrumEnergyRatio;
    % 根据两个时域信号，通过一阶指数加权滑动平均形式计算各自RMS值序列及谱能量比序列，及它们分别的平均值
    Funcs.calcRMSAndSpectrumEnergyRatio_RealTime = @calcRMSAndSpectrumEnergyRatio_RealTime;

    %% 扰动激励生成
    Funcs.generateExcitation = @generateExcitation;

end

%% 获取二自由度DVA的结构参数信息
function [m1,c1,k1,f1,w1,R,Bl,L,m0,c0,k0,f0,w0] = getDVAAndMainSystemParams()
    m1 = 2.69; % 动子质量
    R = 11.197; % 线圈电阻
    L = 2.328e-3; % 线圈电感
    k1 = 31250; % 片状弹簧刚度
    c1 = 50; % 阻尼常数
    Bl = 14;
    w1 = 31.3 * 2 * pi;
    k1 = m1 * w1 ^ 2
    % 求解固有频率
    wn_origin = sqrt(k1 / m1);
    f1 = wn_origin / 2 / pi
    alpha = Bl/R;

    % 主系统参数
    % m0 = 8.259; % 宿主结构质量  -- 不带吸振器
    m0 = 9.607; % 带吸振器附加质量
    c0 = 800;
    % k0 = 1.1965e5;
    % w0 = sqrt(k0 / m0); % 在叠加动子之后主系统的重量会增加，因此固有频率会降低
    % 求解固有频率
    f0 = 30.4
    w0 = f0 * 2 * pi;
    k0 = 4 * pi^2 * f0^2 * m0
end


%% 根据两个时域信号计算其各自RMS值及谱能量比
% 输入：
%   x0 - 主系统响应信号
%   x1 - 吸振器响应信号
% 输出：
%   rmsX0 - x0的RMS值
%   rmsX1 - x1的RMS值
%   spectrumEnergyRatio - Gamma_S = RMS_x1^2 / RMS_x0^2
function [rmsX0, rmsX1, spectrumEnergyRatio] = calcRMSAndSpectrumEnergyRatio(x0, x1)
    rmsX0 = sqrt(mean(x0.^2));
    rmsX1 = sqrt(mean(x1.^2));
    spectrumEnergyRatio = (rmsX1^2) / (rmsX0^2);
end

%% 根据两个时域信号，通过一阶指数加权滑动平均形式计算各自RMS值序列及谱能量比序列，及它们分别的平均值
% 输入：
%   x0 - 主系统响应信号
%   x1 - 吸振器响应信号
%   alpha - 遗忘因子 (0<alpha<=1)，越接近0响应越慢
% 输出：
%   rmsX0Seq - x0的指数加权RMS序列
%   rmsX1Seq - x1的指数加权RMS序列
%   rmsX0Mean - x0的RMS均值
%   rmsX1Mean - x1的RMS均值
%   spectrumEnergyRatioSeq - 指数加权谱能量比序列
%   spectrumEnergyRatioMean - 谱能量比序列的均值
function [rmsX0Seq, rmsX1Seq, rmsX0Mean, rmsX1Mean, spectrumEnergyRatioSeq, spectrumEnergyRatioMean] = calcRMSAndSpectrumEnergyRatio_RealTime(x0, x1, alpha)
    N = length(x0);
    rmsX0Seq = zeros(1, N);
    rmsX1Seq = zeros(1, N);
    spectrumEnergyRatioSeq = zeros(1, N);

    % 初始化
    x0_ewma = x0(1)^2;  % 初值设为第一个样本的平方
    x1_ewma = x1(1)^2;

    % 第一个时刻
    rmsX0Seq(1) = sqrt(x0_ewma);
    rmsX1Seq(1) = sqrt(x1_ewma);
    spectrumEnergyRatioSeq(1) = x1_ewma / x0_ewma;

    % 迭代更新
    for k = 2:N
        % 正确的更新公式
        x0_ewma = (1-alpha)*x0_ewma + alpha*x0(k)^2;
        x1_ewma = (1-alpha)*x1_ewma + alpha*x1(k)^2;

        % 实时RMS
        rmsX0Seq(k) = sqrt(x0_ewma);
        rmsX1Seq(k) = sqrt(x1_ewma);

        % 实时谱能量比
        spectrumEnergyRatioSeq(k) = x1_ewma / x0_ewma;
    end

    % 求均值
    rmsX0Mean = mean(rmsX0Seq);
    rmsX1Mean = mean(rmsX1Seq);
    spectrumEnergyRatioMean = mean(spectrumEnergyRatioSeq);
end

%% 生成不同的扰动激励 (简谐、窄带、随机或白噪声激励)
% 输入：type,freq,bw,fs,totalTime 扰动激励类型，简谐或窄带的中心频率，窄带的相对带宽，采样率，持续时间
% 输出：x 扰动激励信号, t 时间序列
function [x, t] = generateExcitation(type, freq, bw, fs, totalTime, nHarmonics)
    % 生成不同的扰动激励 (简谐、窄带、随机、白噪声或平稳窄带激励)
    % 输入：
    %   type: 'sine', 'narrowband', 'random', 'whitenoise', 'stationary-narrowband'
    %   freq: 中心频率 (Hz)
    %   bw: 相对带宽 (0-1)，窄带激励用
    %   fs: 采样频率 (Hz)
    %   totalTime: 总时间 (s)
    %   nHarmonics: 平稳窄带激励时叠加的简谐波数量（可选）
    % 输出：
    %   x: 激励信号
    %   t: 时间向量
    
    if nargin < 6
        nHarmonics = 200; % 默认叠加200个简谐波
    end

    t = (0:1/fs:totalTime-1/fs);

    switch lower(type)
        case 'sine'
            x = sin(2*pi*freq*t);

        case 'narrowband'
            noise = randn(size(t));
            bwAbs = freq * bw; % 绝对带宽
            fLow = max(0.01, freq - bwAbs/2);
            fHigh = freq + bwAbs/2;
            [b,a] = butter(4, [fLow fHigh]/(fs/2), 'bandpass');
            x = filtfilt(b,a,noise);

        case 'random'
            x = randn(size(t));

        case 'whitenoise'
            x = randn(size(t));

        case 'stationary-narrowband'
            bwAbs = freq * bw; % 绝对带宽
            fLow = max(0.01, freq - bwAbs/2);
            fHigh = freq + bwAbs/2;

            % 随机生成nHarmonics个频率和初相位,使用相同的种子保证结果一致
            rng('default') ;
            rng(42) 
            freqs = fLow + (fHigh - fLow) * rand(1, nHarmonics);
            phases = 2*pi*rand(1, nHarmonics);

            x = zeros(size(t));
            for i = 1:nHarmonics
                x = x + sin(2*pi*freqs(i)*t + phases(i));
            end
            % 归一化幅值
            x = x / sqrt(nHarmonics);

        otherwise
            error('Unknown excitation type.');
    end
end

%% 获取二自由度DVA状态方程的状态矩阵和输入矩阵
function [A,B] = getDVAWithMainSystemStateMatrix(m0, k0, m1, k1, c1)
    A = [0,1,0,0; -(k0+k1)/m0,-(c1)/m0,k1/m0,c1/m0; 0,0,0,1; k1/m1,c1/m1,-k1/m1,-c1/m1];
    B = [0;1/m0;0;0];
end

%% 获取二自由度DVA的f0→a_0,f0→a_1的传递函数(加速度)
function [f02a0,f02a1,a02a1] = calcuTransferAccAmplitude(s, m0, k0, m1, k1, c1)
      % 定义系统的质量、阻尼、刚度矩阵
      M = [m0, 0; 0, m1];
      C = [0, -c1; -c1, c1];
      K = [k0 + k1, -k1; -k1, k1];
  
      % 动态刚度矩阵（频域）
      A = M * s^2 + C * s + K;
  
      % 单位激励作用在主结构（自由度1）上
      F = [1; 0];
  
      % 解频域位移
      X = A \ F;
  
      % 加速度为 s^2 * X
      a0 = abs(s^2 * X(1));
      a1 = abs(s^2 * X(2));
      
      % 输出
      f02a0 = a0;
      f02a1 = a1;
      a02a1 = a1 / a0;
end

%% 获取二自由度DVA的f0→v_0,f0→v_1的传递函数(速度)
function [f02v0,f02v1] = calcuTransferVelocityAmplitude(s, m0, k0, m1, k1, c1)
      % 定义系统的质量、阻尼、刚度矩阵
      M = [m0, 0; 0, m1];
      C = [0, -c1; -c1, c1];
      K = [k0 + k1, -k1; -k1, k1];
  
      % 动态刚度矩阵（频域）
      A = M * s^2 + C * s + K;
  
      % 单位激励作用在主结构（自由度1）上
      F = [1; 0];
  
      % 解频域位移
      X = A \ F;
  
      % 速度为 s * X
      v0 = abs(s * X(1));
      v1 = abs(s * X(2));
      
      % 输出
      f02v0 = v0;
      f02v1 = v1;
end

%% 获取二自由度DVA的x0'→x1'的谱能量比指标(解析解)
% \frac{m_{0} \left(c_{1}^{2} k_{0} + k_{1}^{2} m_{0} + k_{1}^{2} m_{1}\right)}{c_{1}^{2} k_{0} m_{0} + k_{0}^{2} m_{1}^{2} - 2 k_{0} k_{1} m_{0} m_{1} + k_{1}^{2} m_{0}^{2} + k_{1}^{2} m_{0} m_{1}}
function [GammaS] = calcuSpectralEnergyRatioUnderRandomExcite(m0, k0, m1, k1, c1)
    % 计算分子
    numerator = m0 * (c1^2 * k0 + k1^2 * m0 + k1^2 * m1);
    
    % 计算分母
    denominator = c1^2 * k0 * m0 + k0^2 * m1^2 - 2 * k0 * k1 * m0 * m1 + k1^2 * m0^2 + k1^2 * m0 * m1;
    
    % 计算谱能量比指标
    GammaS = numerator / denominator;
end

%% 获取二自由度DVA的f0到x0',x1'的谱能量指标(H2指标)(解析解)
% I_{G0}=\frac{\pi \left(c_{1}^{2} k_{0} m_{0} + k_{0}^{2} m_{1}^{2} - 2 k_{0} k_{1} m_{0} m_{1} + k_{1}^{2} m_{0}^{2} + k_{1}^{2} m_{0} m_{1}\right)}{c_{1} k_{0}^{2} m_{0} m_{1}^{2}}
% I_{G1}=\frac{\pi \left(c_{1}^{2} k_{0} + k_{1}^{2} m_{0} + k_{1}^{2} m_{1}\right)}{c_{1} k_{0}^{2} m_{1}^{2}}
function [IG0, IG1] = calcuSpectralEnergyUnderRandomExcite(m0, k0, m1, k1, c1)
    % 计算分子和分母
    numerator_G0 = c1^2 * k0 * m0 + k0^2 * m1^2 - 2 * k0 * k1 * m0 * m1 + k1^2 * m0^2 + k1^2 * m0 * m1;
    denominator_G0 = c1 * k0^2 * m0 * m1^2;
    
    numerator_G1 = c1^2 * k0 + k1^2 * m0 + k1^2 * m1;
    denominator_G1 = c1 * k0^2 * m1^2;
    
    % 计算谱能量指标
    IG0 = pi * numerator_G0 / denominator_G0;
    IG1 = pi * numerator_G1 / denominator_G1;
end

%% 获取二自由度DVA的x0'→x1'的谱能量比指标关于k1的极大值点(解析解)
% \left(\frac{\omega_{0}^{2}}{2} + \frac{\sqrt{\omega_{0}^{4} + \frac{4 c_{1}^{2} k_{0}}{m_{1}^{2} \left(m_{0} + m_{1}\right)}}}{2}\right)m_1
function [k1_max] = calcuMaximumK1SpectralEnergyRatioUnderRandomExcite(m0, k0, m1, c1)
    % 计算w0
    w0 = sqrt(k0 / m0);
    
    % 计算k1的极大值点
    omega0_squared = w0^2;
    term1 = omega0_squared / 2;
    term2 = sqrt(omega0_squared^2 + (4 * c1^2 * k0) / (m1^2 * (m0 + m1)));
    k1_max_squared = term1 + term2 / 2;
    k1_max = k1_max_squared*m1;
end

%% 获取二自由度DVA的x0→x1的谱能量比指标(解析解)
% \frac{c_{1}^{2} k_{0} m_{0} + c_{1}^{2} k_{0} m_{1} + k_{0} k_{1} m_{1}^{2} + k_{1}^{2} m_{0}^{2} + 2 k_{1}^{2} m_{0} m_{1} + k_{1}^{2} m_{1}^{2}}{c_{1}^{2} k_{0} m_{0} + c_{1}^{2} k_{0} m_{1} + k_{0}^{2} m_{1}^{2} - 2 k_{0} k_{1} m_{0} m_{1} - k_{0} k_{1} m_{1}^{2} + k_{1}^{2} m_{0}^{2} + 2 k_{1}^{2} m_{0} m_{1} + k_{1}^{2} m_{1}^{2}}
function [GammaS] = calcuSpectralEnergyRatioUnderRandomExciteWithDis(m0, k0, m1, k1, c1)
    % 计算分子
    numerator = c1^2 * k0 * m0 + c1^2 * k0 * m1 + k0 * k1 * m1^2 + ...
                k1^2 * m0^2 + 2 * k1^2 * m0 * m1 + k1^2 * m1^2;
    
    % 计算分母
    denominator = c1^2 * k0 * m0 + c1^2 * k0 * m1 + k0^2 * m1^2 - ...
                  2 * k0 * k1 * m0 * m1 - k0 * k1 * m1^2 + ...
                  k1^2 * m0^2 + 2 * k1^2 * m0 * m1 + k1^2 * m1^2;
    
    % 计算谱能量比指标
    GammaS = numerator / denominator;
end

%% 获取二自由度DVA的f0到x0,x1的谱能量指标(H2指标)(解析解)
% I_{G0}=\frac{\pi \left(c_{1}^{2} k_{0} m_{0} + c_{1}^{2} k_{0} m_{1} + k_{0}^{2} m_{1}^{2} - 2 k_{0} k_{1} m_{0} m_{1} - k_{0} k_{1} m_{1}^{2} + k_{1}^{2} m_{0}^{2} + 2 k_{1}^{2} m_{0} m_{1} + k_{1}^{2} m_{1}^{2}\right)}{c_{1} k_{0}^{3} m_{1}^{2}}
% I_{G1}=\frac{\pi \left(c_{1}^{2} k_{0} m_{0} + c_{1}^{2} k_{0} m_{1} + k_{0} k_{1} m_{1}^{2} + k_{1}^{2} m_{0}^{2} + 2 k_{1}^{2} m_{0} m_{1} + k_{1}^{2} m_{1}^{2}\right)}{c_{1} k_{0}^{3} m_{1}^{2}}
function [IG0, IG1] = calcuSpectralEnergyUnderRandomExciteWithDis(m0, k0, m1, k1, c1)
    % 计算分子和分母
    numerator_G0 = c1^2 * k0 * m0 + c1^2 * k0 * m1 + ...
                   k0^2 * m1^2 - 2 * k0 * k1 * m0 * m1 - ...
                   k0 * k1 * m1^2 + k1^2 * m0^2 + ...
                   2 * k1^2 * m0 * m1 + k1^2 * m1^2;
    numerator_G1 = c1^2 * k0 * m0 + c1^2 * k0 * m1 + ...
                    k0 * k1 * m1^2 + k1^2 * m0^2 + ...
                    2 * k1^2 * m0 * m1 + k1^2 * m1^2;
    denominator = c1 * k0^3 * m1^2;

    
    % 计算谱能量指标
    IG0 = pi * numerator_G0 / denominator;
    IG1 = pi * numerator_G1 / denominator;
end

%% 获取二自由度DVA的x0→x1的谱能量比指标关于k1的极大值点(解析解)
function [k1_max] = calcuMaximumK1SpectralEnergyRatioUnderRandomExciteWithDis(m0, k0, m1, c1)
    numA = -2*(m0+m1)^3;
    numB = 2*k0*m1*(m0+m1)^2;
    numC = k0*(2*c1^2*m0^2+4*c1^2*m0*m1+2*c1^2*m1^2+k0*m1^3);

    % 最终的 k1 最大值
    k1_max = (-numB - sqrt(numB^2 - 4*numA*numC))/(2*numA);
end