%function [gammaplus, W, pplus, pminus] = FMM_2D_gamma_e(eps11,eps22,eps33,alpha,beta,h,lambda,N)
function [W, gammaplus, gammaminus] = FMM_2D_gamma_e_new_1_4_matrix(eps11,eps22,eps33,alpha,beta,lambda,N)

k0 = 2*pi/lambda;

%[M] = FMM_2D_matr_eps123(N, alpha, beta, eps11, eps22, eps33, k0);
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
%{
M1 = cat(2, M11, M12, M13, M14);
M2 = cat(2, M21, M22, M23, M24);
M3 = cat(2, M31, M32, M33, M34);
M4 = cat(2, M41, M42, M43, M44);

M = cat(1, M1, M2, M3, M4);
%}
L_EH_up = cat(2,M13,M14);
L_EH_down = cat(2,M23,M24);
L_EH = cat(1,L_EH_up,L_EH_down);

L_HE_up = cat(2,M31,M32);
L_HE_down = cat(2,M41,M42);
L_HE = cat(1,L_HE_up,L_HE_down);

L_full = L_HE*L_EH;
[H_1_4, gamma_sqr_1_4] = eig(L_full);


gamma_1_4 = diag(gamma_sqr_1_4.^0.5);
n_plus_1_4 = 0;
n_minus_1_4 = 0;
NN=(2*N+1)*(2*N+1);
for i=1:2*NN
    if real(gamma_1_4(i))+imag(gamma_1_4(i))>0
        n_plus_1_4 = n_plus_1_4+1;
    else
        n_minus_1_4 = n_minus_1_4 + 1;
        gamma_1_4(i) = -gamma_1_4(i);
    end
end
    
    E_1_4 = L_EH*H_1_4/diag(gamma_1_4);
    %H_1_4 = diag(gamma_1_4)\L_HE*E_1_4;
    W_up = cat(2,E_1_4,E_1_4);
    W_down = cat(2,H_1_4,-H_1_4);
    W = cat(1,W_up,W_down);
    
    gammaplus = gamma_1_4;
    gammaminus = -gamma_1_4;
    
%{
    pplusv = zeros(2*NN,1);
    pminusv = zeros(2*NN,1);
    for m=1:(2*NN)
        pplusv(m) = exp(1i*gammaplus(m)*h);
        pminusv(m) = exp(-1i*gammaminus(m)*h);
    end
    
    pplus = diag(pplusv);
    pminus = diag(pminusv);
%}
