%% 傅里叶变换工具
function Funcs = FuncsFftUtil   
    % % 获取信号幅值谱；可以选择是否去除直流分量及快速展示
    Funcs.GetAmplitudeFreqCurve = @GetAmplitudeFreqCurve;
    % % 获取信号幅值谱及相位谱；可以选择是否去除直流分量及快速展示
    Funcs.GetAmplitudeAndPhaseFreqCurve = @GetAmplitudeAndPhaseFreqCurve;
    % % 低通滤波
    Funcs.LowPassFilter = @LowPassFilter;
    % % 频谱抽取，例如抽取整数附近的频谱
    Funcs.SpectralSampling = @SpectralSampling;
    % % 相位谱处理，所有谱的相位差应该在 ±180°以内
    Funcs.PhaseSpectralProcess = @PhaseSpectralProcess;
    
end

%% 获取信号幅值谱及相位谱；可以选择是否去除直流分量及快速展示
% signal：原始信号
% fs：采样频率
% filtSignal： 滤波后的频率
% removeDcComp: 去掉直流分量,bool值
% fastShow: 是否快速展示
% limit: 幅值谱和相位谱的展示区间
function [Fy,Py] = GetAmplitudeAndPhaseFreqCurve(signal,removeDcComp,fastShow,fs,limit,prefixTitle,figNum)
    limitEnable = true;
    if nargin < 7
        figNum = 0;
    end
    if nargin < 5
        prefixTitle = "";
        limitEnable = false;
    end
    if nargin < 4
        limitEnable = false;
    end
    if nargin < 3
        fs = 1;
    end
    if nargin < 2
        fastShow = false;
    end

    % 幅值谱
    num_of_stemp = length(signal);
    % total_time = len/fs;
    xk = fft(signal); 
    xk_m = abs(xk);
    Fy = xk_m*2/num_of_stemp; %信号的频谱,乘2是因为要忽略负频率分布
    Fy = Fy(1:num_of_stemp/2);
    if removeDcComp
        Fy(1) = 0; %去除直流分量 
    end
    
    % 相位谱
    xkTemp = xk;
    threshold = max(abs(xk))/1000;
    xkTemp(abs(xk)<threshold)=0;
    Py = atan2(imag(xkTemp),real(xkTemp))*180/pi;
    % Py = angle(xk);


    ft = (0:num_of_stemp/2-1)/num_of_stemp*fs;
    T = 1/fs;
    t = 0:T:(num_of_stemp-1)*T;
    if fastShow
        if figNum > 0
            figure(figNum);
        else
            figure();
        end
        clf
        subplot(3,1,1)
        plot(t,signal,'b');
        title(prefixTitle+'时域')
        xlabel('时间（s）')
        ylabel('幅值')
        % xlim([0 total_time])
    
        subplot(3,1,2)
        plot(ft,Fy,'b');
        title(prefixTitle+'幅值谱')
        xlabel('频率（Hz）')
        ylabel('幅值')
        if limitEnable
            xlim(limit);
        end

        subplot(3,1,3)
        Py = Py(1:num_of_stemp/2);
        plot(ft,Py,'b');
        title(prefixTitle+'相位谱')
        xlabel('频率（Hz）')
        ylabel('相位（°）')
        if limitEnable
            xlim(limit);
        end
    end
end

%% 获取信号幅值谱；可以选择是否去除直流分量及快速展示
% signal：原始信号
% fs：采样频率
% filtSignal： 滤波后的频率
% removeDcComp: 去掉直流分量,bool值
% fastShow: 是否快速展示
% limit: 幅值谱和相位谱的展示区间
function [Fy] = GetAmplitudeFreqCurve(signal,removeDcComp,fastShow,fs,limit,prefixTitle)    
    limitEnable = true;
    if nargin < 5
        prefixTitle = "";
        limitEnable = false;
    end
    if nargin < 4
        limitEnable = false;
    end
    if nargin < 3
        fs = 1;
    end
    if nargin < 2
        fastShow = false;
    end

    % 幅值谱
    num_of_stemp = length(signal);
    % total_time = len/fs;
    xk = fft(signal); 
    xk_m = abs(xk);
    Fy = xk_m*2/num_of_stemp; %信号的频谱,乘2是因为要忽略负频率分布
    Fy = Fy(1:floor(num_of_stemp/2));
    if removeDcComp
        Fy(1) = 0; %去除直流分量 
    end

    ft = (0:num_of_stemp/2-1)/num_of_stemp*fs;
    T = 1/fs;
    t = 0:T:(num_of_stemp-1)*T;
    if fastShow
        figure();
        subplot(2,1,1)
        plot(t,signal,'b');
        title(prefixTitle+'时域')
        xlabel('时间（s）')
        ylabel('幅值')
        % xlim([0 total_time])
    
        subplot(2,1,2)
        plot(ft,Fy,'b');
        title(prefixTitle+'幅值谱')
        xlabel('频率（Hz）')
        ylabel('幅值')
        if limitEnable
            xlim(limit);
        end
    end
end

%% 巴特沃斯低通滤波器
% signal：原始信号
% cutOffFreq：截止频率
% fs：采样频率
% filtSignal： 滤波后的频率
function [filtSignal] = LowPassFilter(signal,cutOffFreq,fs)
    order = 4;  % 滤波器阶数
    [b, a] = butter(order, cutOffFreq/(fs/2), 'low');
    filtSignal = filter(b, a, signal);
end

% % 通常的使用方法  -  一般不需要
% fftUtil.GetAmplitudeFreqCurve(fy,false,true,fs);
% fy = fftUtil.LowPassFilter(fy,50,fs);
% figure()
% plot(ft,fy,'b')
% title('驱动信号（滤波后）-幅值谱')
% xlabel('频率（Hz）')
% ylabel('幅值')
% xlim([0,300])
% fftUtil.GetAmplitudeFreqCurve(py,false,true,fs);
% py = fftUtil.LowPassFilter(py,10,fs);
% figure()
% plot(ft,py,'b')
% title('驱动信号-相位图')
% xlabel('频率（Hz）')
% ylabel('相位（°）')
% xlim([0,1000])

%% 频谱抽取&&线性插值，例如抽取整数附近的频谱
% ftTarget：目标频谱X轴，如0:floor(fs/2)
% ft：原始频谱X轴
% fy：原始频谱
% linearInterpolationFt：如果给了这个参数则进行线性插值
function [ftNew,fyNew] = SpectralSampling(ftTarget,ft,fy,linearInterpolationFt)
    % 遍历目标值  
    for i = 1:length(ftTarget)  
        % 计算当前目标值与矩阵中所有元素的绝对差值  
        abs_diffs = abs(ft - ftTarget(i));  
        % 找到最小绝对差值对应的索引  
        [~, idx] = min(abs_diffs);  
        % 将对应的元素存储到结果数组中  
        fyNew(i) = fy(idx);
        ftNew(i) = ft(idx);
    end
    
    if nargin == 4
       y_new = interp1(ftNew, fyNew, linearInterpolationFt, 'linear');
       ftNew = linearInterpolationFt;
       fyNew = y_new;
    end

end

%% 相位谱处理，所有谱的相位差应该某个360°的区间内
function [pyNew] = PhaseSpectralProcess(py,pyMax)
    if nargin < 2
        pyMax = 180;
    end
    pyMin = pyMax - 360;
    for i = 1:length(py)
        pyNew(i) = py(i);
        if py(i) > pyMax
            pyNew(i) = py(i) - 360;
        end
        if py(i)<pyMin
            pyNew(i) = py(i) + 360;
        end
    end
end
