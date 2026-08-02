function host = buildHost(p)
% buildHost  从 getQzsDvaParams 提取宿主结构参数
%
%   host = buildHost(p)
%
%   宿主 SDOF:
%       M*xdd + C*xd + K*x = F(t) + F_react
%
%   字段:
%       M, C, K   - 质量 / 阻尼 / 刚度
%       f0_Hz     - 无阻尼固有频率
%       wn, zeta  - 固有圆频率 / 阻尼比（衰减时间常数 = 1/(zeta*wn)）
%       mu        - 质量比 m_DVA / M

arguments
    p struct
end

host.M     = p.m0;
host.C     = p.c0;
host.K     = p.k0;
host.f0_Hz = p.f0_Hz;
host.wn    = 2 * pi * host.f0_Hz;
host.zeta  = p.xi0;                 % = C/(2*sqrt(M*K))
host.mu    = p.m1 / p.m0;
end
