%function [gammaplus, W, pplus, pminus] = FMM_2D_gamma_e(eps11,eps22,eps33,alpha,beta,h,lambda,N)
function [W, pplus, pminus] = FMM_2D_gamma_e(eps11,eps22,eps33,alpha,beta,h,lambda,N)
%
% function for one layer of periodic structure
% INPUT:
% <epsilon>: vector with permitivity distribution
% <period>: period of the grating in micron
% <lambda>: wavelength in micron
% <theta>: angle of incidence in degree
% <refIndices>: vector [n1 n2] with refr.indices of surrounding
% <N>: number of Fourier orders
%
% OUTPUT:
% <beta>: vector with propagation constants of grating modes
% <e>: electric field vectors
% <W>: interface matrix
% <Pl, Pr>: left and right layer matrix
% program for TE polarization
% for TM polarization we switch eps <-> -mu and E<->H

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

M1 = cat(2, M11, M12, M13, M14);
M2 = cat(2, M21, M22, M23, M24);
M3 = cat(2, M31, M32, M33, M34);
M4 = cat(2, M41, M42, M43, M44);

M = cat(1, M1, M2, M3, M4);



[EH, gamma] = eig(M);
gamma_v = diag(gamma);

NN=(2*N+1)*(2*N+1);


gammaplus = zeros(2*NN,1);
gammaminus = zeros(2*NN,1);
EHplus = zeros(4*NN,2*NN);
EHminus = zeros(4*NN,2*NN);
nplus = 0;
nminus = 0;

for i=1:4*NN
  
   if ( real(gamma(i,i))+imag(gamma(i,i)) )>0
       nplus = nplus + 1;
       gammaplus(nplus) = gamma(i,i);
       EHplus(:,nplus) = EH(:,i);   
   else
       nminus = nminus + 1;
       gammaminus(nminus) = gamma(i,i);
       EHminus(:,nminus) = EH(:,i);
   end        
end
EHp = EHplus(:,1:nplus);
EHm = EHminus(:,1:nminus);
W = cat(2, EHp, EHm);
%gammanew = cat(1,gammaplus,gammaminus);

pplusv = zeros(2*NN,1);
pminusv = zeros(2*NN,1);
for m=1:(2*NN)
    pplusv(m) = exp(1i*gammaplus(m)*h);
    pminusv(m) = exp(1i*gammaminus(m)*h);
end
pplus = diag(pplusv);
pminus = diag(pminusv);


end