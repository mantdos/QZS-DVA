function Funcs = FuncsCommonUtil()
    % 获取一维矩阵最中间的元素
    Funcs.mat_middle_element = @mat_middle_element;
    % 一元二次方程 ax^2+bx+c>0解集计算
    Funcs.equadRootSolve = @equadRootSolve;
    % 计算状态空间方程的龙格库塔方法
    Funcs.rungeKutta = @rungeKutta;
    % 计算带有扰动项状态空间方程的龙格库塔方法
    Funcs.rungeKuttaWithDisturbance = @rungeKuttaWithDisturbance;
    % 数据预处理方法
    Funcs.dataPreprocessing = @dataPreprocessing;
    % 创建一个fig视图，并自动将fig编号增1
    Funcs.createFigure = @createFigure;
    % % 低通滤波
    Funcs.LowPassFilter = @LowPassFilter;
    % % 生成不同类型的激励信号数组
    Funcs.generate_excitation_signal = @generate_excitation_signal;
    
end


function [middle_element] = mat_middle_element(A)
    % 计算矩阵长度
    n = length(A);

    % 计算中间元素的索引
    if mod(n, 2) == 0
        % 如果长度为偶数，取中间两个元素的平均值
        middle_idx = (n/2):(n/2 + 1);
        middle_element = mean(A(middle_idx));
    else
        % 如果长度为奇数，取中间一个元素
        middle_idx = (n + 1) / 2;
        middle_element = A(middle_idx);
    end
end

% 一元二次方程 ax^2+bx+c>0解集计算
function equadRootSolve(a,b,c)
    % 计算判别式
    delta = b^2 - 4*a*c;
    
    % 求解方程（尽管对于不等式来说不是必需的，但可以用于验证）
    roots_of_equation = roots([a, b, c]);
    
    % 确定不等式的解集
    if delta < 0
        disp('不等式的解集是全体实数。');
    elseif delta == 0
        disp(['不等式的解集是除了 x = ', num2str(-b/(2*a)), ' 以外的所有实数。']);
    else
        % 排序根（如果需要）
        roots_of_equation = sort(roots_of_equation);
        disp(['不等式的解集是 x < ', num2str(roots_of_equation(1)), ' 或 x > ', num2str(roots_of_equation(2)), '。']);
    end
end

% 用于计算状态空间方程的龙格库塔方法
function [x,xDot] = rungeKutta(h,A,B,x0,u0)
    % 计算第一阶导数
    k1 = A*x0 + B*u0;

    % 计算第二阶导数
    k2 = A*(x0 + h/2*k1)+B*u0;        

    % 计算第三阶导数
    k3 = A*(x0 + h/2*k2)+B*u0;          

    % 计算第四阶导数    
    k4 = A*(x0 + h*k3)+B*u0;              

    % 计算状态空间方程的解
    x = x0 + h/6*(k1 + 2*k2 + 2*k3 + k4);
    xDot = k1;

    % % 计算第一阶导数
    % k1 = A*x0 + B*u0;
    % 
    % % 计算第二阶导数
    % k2 = A*(x0 + h/4*k1)+B*u0;        
    % 
    % % 计算第三阶导数
    % k3 = A*(x0 + h/8*k1 + h/8*k2)+B*u0;          
    % 
    % % 计算第四阶导数    
    % k4 = A*(x0 - h/2*k2 + h*k3)+B*u0;      
    % 
    % % 计算第五阶导数    
    % k5 = A*(x0 + 3*h/16*k1 + 9*h/16*k4)+B*u0; 
    % 
    % % 计算第六阶导数    
    % k6 = A*(x0 - 3*h/7*k1 + 2*h/7*k2 + 12*h/7*k3 -12*h/7*k4 + 8*h/7*k5)+B*u0;          
    % 
    % % 计算状态空间方程的解
    % x = x0 + h/90*(7*k1 + 32*k3 + 12*k4 + 32*k5 + 7*k6);
    % xDot = k1;
end

% 计算带有扰动项状态空间方程的龙格库塔方法
function [x,xDot] = rungeKuttaWithDisturbance(h,A,B,H,x0,u0,d0)
    % 计算第一阶导数
    k1 = A*x0 + B*u0 + H*d0;

    % 计算第二阶导数
    k2 = A*(x0 + h/2*k1)+B*u0 + H*d0;      

    % 计算第三阶导数
    k3 = A*(x0 + h/2*k2)+B*u0 + H*d0;         

    % 计算第四阶导数    
    k4 = A*(x0 + h*k3)+B*u0 + H*d0;            

    % 计算状态空间方程的解
    x = x0 + h/6*(k1 + 2*k2 + 2*k3 + k4);  
    xDot = k1;
    
    % % 计算第一阶导数
    % k1 = A*x0 + B*u0 + H*d0;
    % 
    % % 计算第二阶导数
    % k2 = A*(x0 + h/4*k1)+B*u0 + H*d0;      
    % 
    % % 计算第三阶导数
    % k3 = A*(x0 + h/8*k1 + h/8*k2)+B*u0 + H*d0;        
    % 
    % % 计算第四阶导数    
    % k4 = A*(x0 - h/2*k2 + h*k3)+B*u0 + H*d0;  
    % 
    % % 计算第五阶导数    
    % k5 = A*(x0 + 3*h/16*k1 + 9*h/16*k4)+B*u0 + H*d0;
    % 
    % % 计算第六阶导数    
    % k6 = A*(x0 - 3*h/7*k1 + 2*h/7*k2 + 12*h/7*k3 -12*h/7*k4 + 8*h/7*k5)+B*u0 + H*d0;         
    % 
    % % 计算状态空间方程的解
    % x = x0 + h/90*(7*k1 + 32*k3 + 12*k4 + 32*k5 + 7*k6);
    % xDot = k1;
end

% 仿真结果数据预处理
function [s,ft,fs] = dataPreprocessing(signal,total_time,cutTime,t_start,t_end,needRemoveDc)
    num_of_stemp = length(signal);
    T = total_time/(num_of_stemp-1);
    fs = 1/T;
    ft = (0:num_of_stemp/2-1)/num_of_stemp*fs;
    t = 0:T:(num_of_stemp-1)*T;

    if cutTime
        indexStart = floor(t_start*fs);
        indexEnd = floor(t_end*fs);
        total_time = t_end - t_start;
        signal = signal(indexStart:1:indexEnd);
        num_of_stemp = length(signal);
        ft = (0:num_of_stemp/2-1)/num_of_stemp*fs;
        t = 0:T:(num_of_stemp-1)*T;
    end

    if needRemoveDc
        signal = signal - sum(signal)/length(signal);
    end
    
    s = signal;
end

function figN = createFigure(figN)
    figure(figN);
    figN = figN + 1;
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


%% generate_excitation_signal - 生成不同类型的激励信号数组
%
% 语法: excitationArr = generate_excitation_signal(type, params)
%
% 输入:
%   type - 字符串，指定信号类型。可选值:
%          'sine'      : 单频正弦信号
%          'random'    : 宽带随机信号 (白噪声)
%          'chirp'     : 线性扫频信号 (Chirp)
%          'narrow'    : 窄带随机信号
%
%   params - 结构体，包含生成信号所需的参数:
%       .L          - (必需) 信号长度 (点数)
%       .Ts         - (必需) 采样周期 (s)
%       .amplitude  - (必需) 信号幅值
%
%       对于 'sine':
%       .f_sine     - (必需) 正弦信号频率 (Hz)
%
%       对于 'chirp':
%       .f_start    - (必需) 起始频率 (Hz)
%       .f_end      - (必需) 终止频率 (Hz)
%       .T_sweep    - (必需) 扫频持续的总时间 (s)
%
%       对于 'narrow':
%       .f_low      - (必需) 窄带通带的下限频率 (Hz)
%       .f_high     - (必需) 窄带通带的上限频率 (Hz)
%
% 输出:
%   excitationArr - 1xL 的行向量，生成的激励信号
function excitationArr = generate_excitation_signal(type, params)
    % 检查必需的参数
    if ~isfield(params, 'L') || ~isfield(params, 'Ts') || ~isfield(params, 'amplitude')
        error('参数 "params" 必须包含字段: L, Ts, amplitude');
    end

    % 生成时间向量
    t = (0:params.L-1) * params.Ts;
    fs = 1/params.Ts;

    % 使用 switch 根据类型生成信号
    switch lower(type)
        case 'sine'
            % --- 单频正弦信号 ---
            if ~isfield(params, 'f_sine')
                error('对于 "sine" 类型, params 必须包含字段: f_sine');
            end
            f_sine = params.f_sine;
            excitationArr = params.amplitude * sin(2 * pi * f_sine * t);
            
            plot_signal(excitationArr, fs, sprintf('Sine Wave (%.1f Hz)',f_sine));
        case 'random'
            % --- 宽带随机信号 (零均值白噪声) ---
            % 使用 randn 生成高斯白噪声, 效果通常比 rand 好
            rng(100); % 您可以将0替换为任何其他整数来设置不同的随机种子
            noise = randn(1, params.L);
            % 标准化为单位方差，然后乘以幅值
            excitationArr = params.amplitude * (noise / std(noise));

            plot_signal(excitationArr, fs, 'Broadband Random');

        case 'chirp'
            % --- 线性扫频信号 ---
            if ~isfield(params, 'f_start') || ~isfield(params, 'f_end') || ~isfield(params, 'T_sweep')
                error('对于 "chirp" 类型, params 必须包含字段: f_start, f_end, T_sweep');
            end
            % 使用 MATLAB 内置的 chirp 函数
            excitationArr = params.amplitude * chirp(t, params.f_start, params.T_sweep, params.f_end, 'linear');

            plot_signal(excitationArr, fs, sprintf('Chirp Signal (%.1f-%.1f Hz)',params.f_start,params.f_end));

        case 'narrow'
            % --- 窄带随机信号 ---
            if ~isfield(params, 'f_low') || ~isfield(params, 'f_high')
                error('对于 "narrow" 类型, params 必须包含字段: f_low, f_high');
            end
            
            % 1. 生成宽带高斯白噪声
            noise = randn(1, params.L);
            
            % 2. 设计一个带通滤波器 (这里使用一个4阶巴特沃斯滤波器)
            filterOrder = 4;
            bpFilt = designfilt('bandpassiir', 'FilterOrder', filterOrder, ...
                            'HalfPowerFrequency1', params.f_low, ...
                            'HalfPowerFrequency2', params.f_high, ...
                            'SampleRate', fs);
                        
            % 3. 对噪声进行滤波
            filtered_noise = filter(bpFilt, noise);
            
            % 4. 标准化幅值
            % 滤波后信号的能量会改变, 重新标准化使其幅值可控
            excitationArr = params.amplitude * (filtered_noise / max(abs(filtered_noise)));

            plot_signal(excitationArr, fs, sprintf('Narrow-band Random (%.1f-%.1f Hz)',params.f_low,params.f_high));

        otherwise
            error('未知的激励信号类型: %s. 可选值为 ''sine'', ''random'', ''chirp'', ''narrow''', type);
    end
end

%% --- 辅助绘图函数 ---
function plot_signal(signal, fs, title_text)
    L = length(signal);
    t = (0:L-1) / fs;
    
    figure(100);
    clf
    
    % 绘制时域图像
    subplot(2,1,1);
    plot(t, signal);
    title([title_text ' - Time Domain']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    
    % 绘制频域图像 (FFT)
    subplot(2,1,2);
    Y = fft(signal);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(L/2))/L;
    plot(f, P1);
    title([title_text ' - Frequency Domain (FFT)']);
    xlabel('Frequency (Hz)');
    ylabel('|P1(f)|');
    grid on;
    xlim([0, fs/2]); % 只显示到奈奎斯特频率
end