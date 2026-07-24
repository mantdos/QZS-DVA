function plant = buildPlant(p)
% buildPlant  从 getQzsDvaParams 映射到 RFE-AIC 对象参数

arguments
    p struct
end

plant.m  = p.m1;
plant.c  = p.c1;
plant.k1 = p.kn1;          % 残余线性刚度
plant.k3 = p.kn3;          % 三次非线性
plant.k5 = 0;              % 第一版无五次项（预留）
plant.k_m = p.kn1;         % 参考模型刚度 = 残余线性刚度
plant.c_m = p.c1;          % 参考模型阻尼
plant.f1_Hz = p.f1Residual_Hz;
plant.m0 = p.m0;
plant.k0 = p.k0;
plant.c0 = p.c0;
end
