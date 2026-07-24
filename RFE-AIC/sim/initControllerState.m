function st = initControllerState(fltSys, ctl)
% initControllerState  初始化控制器滤波器状态与自适应参数

arguments
    fltSys struct
    ctl    struct
end

nxBp = fltSys.bp.nx;
nxLp = fltSys.lp.nx;
zBp = zeros(nxBp, 1);
zLp = zeros(nxLp, 1);

chans = {'e','ed','edd','y','yd'};
for i = 1:numel(chans)
    st.(['bp_' chans{i}]) = zBp;
    st.(['lp_' chans{i}]) = zLp;
end

if ctl.adapt_on
    st.dk_hat = 0;
    st.dc_hat = 0;
else
    st.dk_hat = ctl.dk_frozen;
    st.dc_hat = ctl.dc_frozen;
end
end
