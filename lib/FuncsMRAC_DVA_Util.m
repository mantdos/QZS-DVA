function Funcs = FuncsMRAC_DVA_Util()
    % 获取二自由度DVA的f0→x_0传递函数
    Funcs.calcuTransferAmplitude = @calcuTransferAmplitude;
    % 获取单独的DVA状态方程的状态矩阵，输入矩阵，扰动输入矩阵
    Funcs.getDVASystemStateMatrix = @getDVASystemStateMatrix;
    % 获取单独的DVA的仿真结果
    Funcs.DVASystemRuggeKuttaSimu = @DVASystemRuggeKuttaSimu;
    % 获取二自由度DVA状态方程的状态矩阵，输入矩阵，扰动输入矩阵
    Funcs.getDVAMainSystemStateMatrix = @getDVAMainSystemStateMatrix;
    % 获取二自由度DVA系统的仿真结果
    Funcs.DVAMainSystemRuggeKuttaSimu = @DVAMainSystemRuggeKuttaSimu;
end

%% 获取二自由度DVA的f0→x_0,f0→y,x_0→y传递函数，传递函数矩阵元素为x0,y(主结构位移，吸振器动质量相对于主结构的位移) 无控制时
function [f02x0,f02y,x02y] = calcuTransferAmplitude(s, m0, k0, m1, k1, c1)
    mat_00 = m0 * s ^ 2 + k0;
    mat_01 = -(c1) * s - k1;
    mat_11 = m1 * s ^ 2;
    mat_12 = m1 * s ^ 2 + (c1) * s + k1;
    mat = [mat_00, mat_01; mat_11, mat_12];
    % f02x0 = abs([1, 0, 0] * s ^ 2 * inv(mat) * [1; 0; 0]);
    f02x0 = abs([1,0]*inv(mat)*[1;0;]);
    f02y = abs([0,1]*inv(mat)*[1;0;]);
    x02y = f02y/f02x0;
end

% 获取单独的DVA状态方程的状态矩阵，输入矩阵，扰动输入矩阵
function [A,B,H] = getDVASystemStateMatrix(m1, k1, c1, Bl, R)
    A = [0, 1; -k1/m1, -c1/m1];
    B = [0; 1/m1*Bl/R]; % 输入矩阵
    H = [0; -1]; % 扰动输入矩阵
end

% 获取单独的DVA的仿真结果
function [XNext,XDotNext] = DVASystemRuggeKuttaSimu(h, m1, k1, c1, X, u0, d0, Bl, R)
    [A,B,H] = getDVASystemStateMatrix(m1, k1, c1, Bl, R);

    % 计算第一阶导数
    k1 = A*X + B*u0 + H*d0;
        
    % 计算第二阶导数
    k2 = A*(X + h/2*k1)+B*u0 + H*d0;      
    
    % 计算第三阶导数
    k3 = A*(X + h/2*k2)+B*u0 + H*d0;         
    
    % 计算第四阶导数    
    k4 = A*(X + h*k3)+B*u0 + H*d0;            
    
    % 计算状态空间方程的解
    XNext = X + h/6*(k1 + 2*k2 + 2*k3 + k4);  
    XDotNext = k1;
end

% 获取二自由度DVA状态方程的状态矩阵，输入矩阵，扰动输入矩阵
function [A,B,H] = getDVAMainSystemStateMatrix(m0, k0, m1, k1, c1, Bl, R)
    A = [0,1,0,0; -k0/m0,0,k1/m0,c1/m0; 0,0,0,1; k0/m0,0,-(m1*k1+m0*k1)/(m0*m1),-(m1*c1+m0*c1)/(m0*m1)];
    B = [0;-1/m0*Bl/R;0;(m0+m1)/(m0*m1)*Bl/R];
    H = [0;1/m0;0;-1/m0];
end

% 获取二自由度DVA的仿真结果
% X：状态变量，u0：控制力，d0：激励力/扰动
% XNext:[x0,x0_dot,y,y_dot],XDotNext:[x0_dot,x0_ddot,y_dot,y_ddot]
function [XNext,XDotNext] = DVAMainSystemRuggeKuttaSimu(h, m0, k0, m1, k1, c1, X, u0, d0, Bl, R)
    [A,B,H] = getDVAMainSystemStateMatrix(m0, k0, m1, k1, c1, Bl, R);

    % 计算第一阶导数
    k1 = A*X + B*u0 + H*d0;
        
    % 计算第二阶导数
    k2 = A*(X + h/2*k1)+B*u0 + H*d0;      
    
    % 计算第三阶导数
    k3 = A*(X + h/2*k2)+B*u0 + H*d0;         
    
    % 计算第四阶导数    
    k4 = A*(X + h*k3)+B*u0 + H*d0;            
    
    % 计算状态空间方程的解
    XNext = X + h/6*(k1 + 2*k2 + 2*k3 + k4);  
    XDotNext = k1;
end