function [gamma, W, pplus, pminus, eps] = FMM_1D_TE_RT_beta_e(alpha, beta, epsilon, periodx, periody, h, lambda, theta, phi, refIndices, N, M)
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

eps = FMM_happy_epsilon(N,epsilon);
%{
[M, G] = FMM_2D_matr(N, alpha, beta, eps, k0);
[E, gamma2] = eig(M);
gamma = sqrt(gamma2);
H = (1/k0)*gamma\G*E;
Econd = rcond(E);


W_up = cat(2, E, E);
W_down = cat(2, H, -H);
W = cat(1, W_up, W_down);
Wcond = rcond(W);

NN = (2*N+1)*(2*N+1);
pplusv = zeros(2*NN,1);
pminusv = zeros(2*NN,1);
nplus=0;

for m=1:(2*NN)
    pplusv(m) = exp(1i*gamma(m,m)*h);
    pminusv(m) = exp(-1i*gamma(m,m)*h);
end
pplus = diag(pplusv);
pminus = diag(pminusv);

for i=1:(2*NN)
    if (real(gamma(i,i))+imag(gamma(i,i)))>0
        nplus = nplus+1;
    end
end
nplusg = nplus
%}
[M] = FMM_2D_matr_new(N, alpha, beta, eps, k0);
[EH, gamma] = eig(M);
gamma_v = diag(gamma);

NN=(2*N+1)*(2*N+1);
%{
for i=1:(4*NN-1)          
    for j = 1:(4*NN-i) 
      if imag(gamma_v(j))<imag(gamma_v(j+1))%( a[j-1] > a[j] )
          
      x = gamma_v(j);         %x=a[j-1];
      gamma_v(j) = gamma_v(j+1);%a[j-1]=a[j];
      gamma_v(j+1) = x;           %a[j]=x;
      
      xx = EH(:,j);
      EH(:,j) = EH(:,j+1);
      EH(:,j+1) = xx;
      
      end
    end
end
%}

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


gammatot = cat(1,gammaplus,gammaminus);

pplusv = zeros(2*NN,1);
pminusv = zeros(2*NN,1);
for m=1:(2*NN)
    pplusv(m) = exp(1i*gammaplus(m)*h);
    pminusv(m) = exp(1i*gammaminus(m)*h);
end
pplus = diag(pplusv);
pminus = diag(pminusv);

%{
mzero = zeros((2*N+1)*(2*N+1), (2*N+1)*(2*N+1));
miden = eye((2*N+1)*(2*N+1), (2*N+1)*(2*N+1));

Pl_up = cat(2, miden, mzero);
Pl_down = cat(2, mzero, pplus);
Pl = cat(1, Pl_up, Pl_down);

Pr_up = cat(2,pplus, mzero);
Pr_down = cat(2, mzero, miden);
Pr = cat(1, Pr_up, Pr_down);

%p_up = cat(2, pplus, mzero);
%p_down = cat(2, mzero, pminus);
%P = cat(1, p_up, p_down);

%eigenvalues gamma - kz in periodic layer
%eigenvectors e - electric field amplitudes E
%}

end