function st = initControllerState(fltSys, ctl)
% initControllerState  初始化控制器滤波器状态与自适应参数

arguments
    fltSys struct
    ctl    struct
end

nxBp = fltSys.bp.nx; % 带通滤波器的状态变量数量 2
nxLp = fltSys.lp.nx; % 低通滤波器的状态变量数量 1
zBp = zeros(nxBp, 1); % 带通滤波器的状态变量初始值 0
zLp = zeros(nxLp, 1); % 低通滤波器的状态变量初始值 0

% 初始化带通滤波器状态，只过 e、ed、y、yd
chans = {'e','ed','y','yd'};
for i = 1:numel(chans)
    st.(['bp_' chans{i}]) = zBp;
end

% 初始化带通+低通滤波器状态，只过 edd
chans = {'edd'};
for i = 1:numel(chans)
    st.(['bp_' chans{i}]) = zBp;
    st.(['lp_' chans{i}]) = zLp;
end

% 初始化自适应参数
if ctl.adapt_on
    st.dk_hat = 0;
    st.dc_hat = 0;
else
    st.dk_hat = ctl.dk_frozen;
    st.dc_hat = ctl.dc_frozen;
end
end
