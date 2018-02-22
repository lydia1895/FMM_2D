function [M] = FMM_2D_matr_eps123(N, alpha, beta, eps11, eps22, eps33, k0)

NN = (2*N+1)*(2*N+1);
miden = eye(NN,NN);
mzero = zeros(NN,NN);
M11 = mzero;
M12 = mzero;
M13 = (1/k0)*alpha/eps33*beta;
M14 = k0*miden - (1/k0)*alpha/eps33*alpha;

M21 = mzero;
M22 = mzero;
M23 = -k0*miden + (1/k0)*beta/eps33*beta;
M24 = -(1/k0)*beta/eps33*alpha;

M31 = -(1/k0)*alpha*beta;
M32 = -k0*eps22 + (1/k0)*alpha*alpha;
M33 = mzero;
M34 = mzero;

M41 = k0*eps11 -(1/k0)*beta*beta;
M42 = (1/k0)*beta*alpha;
M43 = mzero;
M44 = mzero;

M1 = cat(2, M11, M12, M13, M14);
M2 = cat(2, M21, M22, M23, M24);
M3 = cat(2, M31, M32, M33, M34);
M4 = cat(2, M41, M42, M43, M44);

M = cat(1, M1, M2, M3, M4);