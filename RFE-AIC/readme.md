RFE-AIC： Reaction-Force-Error-driven Adaptive Impedance Control
### 中文名称： 基于反作用力误差驱动的自适应阻抗控制


### ODE仿真代码步骤
#### 1、QZS仿真参数
    需要将../QZS-DVA/getQzsDvaParams.m摞到../lib中作为公共函数获取QZS的所有参数
#### 2、定义仿真所需额外参数
* 采样率 fs
* 算法的控制参数 c_1、c_2、Gamma_k、Gamma_c、omega_c
* 频带选择参数：带通滤波器的相关参数（二阶）
* 频带选择参数、omega_c如何选择，权衡？
#### 3、定义激励
* 随机激励打底的窄带主导激励
* 相关参数：窄带带宽、主导频率、随机激励相对窄带的占比、整体信号的幅值
* 窄带带宽为0、随机激励占比为0时表示简谐激励
#### 4、控制算法运行,注意以下事项：
* \ddot{x}、e、\dot{e}、\ddot{e}都要经过频带选择，然后\ddot{e}额外在经过一阶低通滤波
* 仿真方法采用ODE定步长，还是自己用4阶龙格库塔算法？


# QZS-DVA 自适应阻抗控制仿真（第一版）
## 0. 总体说明

实现一个准零刚度动力吸振器（QZS-DVA）在窄带随机主导激励下的自适应阻抗控制仿真。
第一版目标：验证方法可行性，**不实现投影算子和 σ-modification**，使用裸自适应律，
以便观察参数漂移等原始行为。

代码组织为 MATLAB 函数化脚本，模块清晰、参数集中，便于后续扫参与扩展。

## 1. 目录与公共函数
QZS-DVA/
lib/
    getQzsDvaParams.m      % 从 ../QZS-DVA-basic/ 迁移过来
RFE-AIC/
  funcs/
    genExcitation.m
    envelopeStats.m
    designBandpass.m
    equivStiffnessTheory.m
    rk4Step.m
  sim/
    simDvaOnly.m
    controllerUpdate.m
    runFrozenSweep.m
    runAdaptive.m
  main.m
    main.m
有需要打包出来的公共函数就放在这个文件夹下

## 2. 参数配置（集中在 main.m 顶部结构体）
%% 物理与参考模型
p    = getQzsDvaParams();
。。。

%% 时间与采样
sim.fs_ctrl = 2000;              % 控制器采样率 Hz
sim.M       = 8;                 % 每控制周期内对象积分步数
sim.h       = 1/(sim.fs_ctrl*sim.M);   % 对象积分步长
sim.T       = 120;               % 仿真时长 s
sim.t_discard = 10;              % 前 10 s 不计入指标（排除滤波器瞬态）

%% 频带选择（二阶带通，相对阶必须为 1）
flt.omega_d  = 2*pi*2.0;         % 目标主导频率 rad/s
flt.zeta_b   = 0.1;              % 带通阻尼比，带宽 = 2*zeta_b*omega_d
flt.omega_c  = 10*flt.omega_d;   % 一阶低通截止频率

%% 控制参数
ctl.c1     = 20;                 % 状态误差面系数
ctl.lambda = 1.0;                % 反作用力误差权重（c2 = lambda/omega_c）
ctl.c2     = ctl.lambda/flt.omega_c;
ctl.Ms     = 1 + ctl.c2*flt.omega_c;   % = 1 + lambda
ctl.K      = 50;                 % 趋近增益
ctl.Gamma_k = 1e3;               % 等效刚度自适应增益
ctl.Gamma_c = 1e2;               % 等效阻尼自适应增益
ctl.adapt_on = true;             % false = 参数冻结模式
ctl.dk_frozen = 0;               % 冻结模式下的 Δk_eq
ctl.dc_frozen = 0;               % 冻结模式下的 Δc_eq
ctl.u_sat  = Inf;                % 第一版不加饱和

%% 激励
exc.f_d        = 2.0;            % 主导频率 Hz
exc.bw         = 0.4;            % 窄带带宽 Hz（=0 时为纯简谐）
exc.rand_ratio = 0.3;            % 随机分量功率占比（=0 时为纯简谐）
exc.amp        = 1.0;            % 整体幅值标度
exc.seed       = 1;

**参数化说明**：`c2 = lambda/omega_c` 使 `Ms = 1+lambda` 与 `omega_c` 解耦，
扫参时只扫 `lambda`，避免 `omega_c` 变化间接改变反作用力误差项权重。

## 3. 激励生成（RFE-AIC/funcs/genExcitation.m）

输入：exc 结构体、时间向量 t
输出：X（宿主结构的扰动力）

实现：
1. 主导正弦分量：`A_h*sin(2*pi*f_d*t)`
2. 随机分量：白噪声经中心频率 f_d、带宽 bw 的带通滤波
3. 按 `rand_ratio` 分配两者功率，总功率由 `amp` 标定
4. `bw==0 && rand_ratio==0` 时退化为纯简谐
5. 使用 `rng(exc.seed)` 保证可复现

**同时输出激励的实测统计**：均方值、功率谱（用于校核）。

## 4. 包络统计诊断（RFE-AIC/funcs/envelopeStats.m）

输入：信号 y、时间向量 t、丢弃时长
输出结构体：
- `A_mean`：包络均值
- `A_cv`：变异系数 σ_A/Ā
- `tau_env`：包络自相关衰减到 1/e 的时间
- `f_env_bw`：包络功率谱的等效带宽

实现：Hilbert 变换取模得包络 → 统计量 → 自相关/谱分析。

**每次仿真都调用并打印结果**，用于判断：
- A_cv < 10%：主导正弦占优，等效刚度近似定值
- A_cv > 30%：随机显著，存在固有瞬时失配
- tau_env 与自适应时间常数、带通时间常数（≈2/(2*zeta_b*omega_d)）是否同量级

## 5. 滤波器设计（RFE-AIC/funcs/designBandpass.m）

设计**二阶带通**（相对阶必须为 1）：

```
H_bp(s) = 2*zeta_b*omega_d*s / (s^2 + 2*zeta_b*omega_d*s + omega_d^2)
```

用 Tustin 双线性变换（预畸变到 omega_d）离散化到 fs_ctrl，返回状态空间形式 (Ad,Bd,Cd,Dd)。

**三路同构要求**：e、ė、ë 各过一份**结构完全相同**的带通（三个独立实例，
各自维护状态）。ë 那一路带通之后再串一阶低通 ω_c。

**为保持三路相位严格一致，e 和 ė 两路也串同样的一阶低通。**
即三路都是：带通 → 一阶低通。

输出时打印在 omega_d 处的总相位滞后，应 < 20°。

## 6. 被控对象积分（RFE-AIC/funcs/rk4Step.m + RFE-AIC/sim/simDvaOnly.m）

对象方程（DVA 单独，xdd 为给定输入）：
```
m*ydd + c*yd + k1*y + k3*y^3 + k5*y^5 = -m*xdd + u
```

参考模型：
```
m*ym_dd + c_m*ym_d + k_m*ym = -m*xdd
```

**双时钟结构**：
```matlab
for k = 1:N_ctrl
    u_k = controllerUpdate(...);          % 离散更新，零阶保持
    for i = 1:M
        state = rk4Step(state, u_k, xdd_interp, h);
    end
    % 记录数据
end
```

对象状态：[y; yd; ym; ym_d]，共 4 维。

## 7. 控制器更新（RFE-AIC/sim/controllerUpdate.m）

每个控制周期执行一次，内部维护滤波器状态与自适应参数。

### 7.1 误差量
```
e   = y - ym
ed  = yd - ym_d
edd = ydd - ym_dd     % 由对象方程和参考模型方程直接算，非数值微分
```

### 7.2 带内提取（三路同构）
```
Be   = lowpass(bandpass(e))
Bed  = lowpass(bandpass(ed))
a_f  = lowpass(bandpass(edd))
```
同时对 y、yd 做相同处理得到 `By`、`Byd`（用于回归向量）。

### 7.3 复合误差面
```
s = Bed + c1*Be + c2*a_f
```

### 7.4 已知补偿力
```
F_known = -c_m*Bed - k_m*Be
```

### 7.5 未知集总力估计
```
F_unc_hat = dk_hat*By + dc_hat*Byd
```

### 7.6 控制律
```
u = (m/Ms) * ( -K*s + c2*omega_c*a_f - c1*Bed ) - F_known + F_unc_hat
```
（推导见文档 2.5.5 节，注意符号一致性，实现时逐项核对）

### 7.7 自适应律（第一版：裸形式，无投影、无 σ、无死区）
```
dk_hat = dk_hat - Gamma_k * s * By  * Ts
dc_hat = dc_hat - Gamma_c * s * Byd * Ts
```
使用前向欧拉离散。

**冻结模式**：`ctl.adapt_on == false` 时，`dk_hat = ctl.dk_frozen`，
`dc_hat = ctl.dc_frozen`，不更新。

### 7.8 需要记录的量
每步保存：t, y, yd, ydd, ym, e, ed, edd, s, a_f, u, dk_hat, dc_hat, By, Byd

## 8. 等效刚度理论值（RFE-AIC/funcs/equivStiffnessTheory.m）

第一版实现两个版本供对比：

**版本 A — 单谐波描述函数**（假设纯简谐响应，幅值 A）：
```
dk_harmonic = (3/4)*k3*A^2 + (5/8)*k5*A^4
```

**版本 B — 高斯等效线性化**（假设零均值高斯响应，方差 σ_y²）：
```
dk_gaussian = 3*k3*sigma_y^2 + 15*k5*sigma_y^4
```

两者都返回，供后续对比。

**注意**：混合激励（主导正弦+随机）下严格系数介于两者之间，
依赖确定/随机分量功率比。第一版先用这两个版本框住范围，
后续再实现混合情形的严格公式。

计算时 A 与 σ_y 均由**实际仿真响应**统计得到（丢弃前 t_discard 秒）。

## 9. 实验脚本

### 9.1 runFrozenSweep.m —— 参数冻结天花板测试
**这是第一个要跑的实验。**

流程：
1. 扫描激励幅值 `amp`，跨 1~2 个量级（如 logspace，10 个点）
2. 每个幅值下，先跑一次**无控**仿真（u=0），统计响应均方 σ_y²
3. 用 equivStiffnessTheory 算出 dk_frozen（A、B 两个版本各算一次）
4. 用冻结参数跑受控仿真
5. 每个工况跑 10 个随机种子，结果取统计

输出图：
- **图 1（核心）**：横轴激励幅值（对数），纵轴反共振频率。三条线：无控、
  冻结版本 A、冻结版本 B。**无控应明显上翘，冻结受控应尽量平**
- 图 2：横轴激励幅值，纵轴宿主传递率峰值 / 响应 RMS
- 图 3：每个幅值下的包络统计 A_cv

**判据**：若冻结最优也压不平图 1，说明等效阻抗集总假设不足，
需先改结构再谈自适应。

### 9.2 runAdaptive.m —— 自适应性能测试
天花板确认后再跑。

场景一：固定幅值，观察参数收敛
- 输出：dk_hat、dc_hat 时程（**注意观察是否漂移**）、s 时程、e 时程
- 与冻结最优值对比：dk_hat 是否收敛到该值附近

场景二：幅值阶跃（如 t=60s 时 ×3）
- 输出：反共振频率偏移量、恢复时间、宿主响应超调
- 对比：无控 / 冻结（按初始幅值设定，故阶跃后失配）/ 自适应

场景三：lambda 扫描
- lambda 从 0（纯状态跟踪）扫到较大值
- 输出：性能指标 vs lambda，**检查是否存在中间最优**
- 这张图决定反作用力误差项是否真有价值

### 9.3 指标计算
统一在一个函数里实现，所有实验共用：
- 反共振频率：对受控系统做扫频或用响应谱峰识别
- 宿主响应 RMS 衰减量
- ‖e‖_rms、‖ë‖_rms（**两者都要报**，用于验证"状态收敛 ≠ 加速度收敛"）
- 控制力 u 的 RMS，及其与被动恢复力 RMS 的比值

## 11. 执行顺序

1. 先跑 genExcitation + envelopeStats，检查 A_cv 与 tau_env，**记录结果**
2. 跑无控仿真，复现幅值-反共振频率漂移曲线（对应 1.93→2.5 Hz）
3. 跑 runFrozenSweep，做天花板测试
4. 天花板通过后，跑 runAdaptive 三个场景

每一步的结果保存为 .mat，图保存为 .fig 和 .png。