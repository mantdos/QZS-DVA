function out = equivStiffnessTheory(k3, k5, A, sigmaY)
% equivStiffnessTheory  等效刚度增量理论值（谐波 / 高斯两版）
%
%   out = equivStiffnessTheory(k3, k5, A, sigmaY)
%
%   版本 A — 单谐波描述函数（幅值 A）:
%       dk_harmonic = (3/4)*k3*A^2 + (5/8)*k5*A^4
%   版本 B — 高斯等效线性化（方差 sigmaY^2）:
%       dk_gaussian = 3*k3*sigmaY^2 + 15*k5*sigmaY^4
%
%   混合激励下严格值介于两者之间；第一版用二者框定范围。

arguments
    k3     (1,1) double
    k5     (1,1) double
    A      (1,1) double
    sigmaY (1,1) double
end

out.A = A;
out.sigmaY = sigmaY;
out.sigmaY2 = sigmaY^2;
out.dk_harmonic = (3/4)*k3*A^2 + (5/8)*k5*A^4;
out.dk_gaussian = 3*k3*sigmaY^2 + 15*k5*sigmaY^4;
end
