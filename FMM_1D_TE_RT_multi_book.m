function [gamma, eta_R, eta_T, eps] = FMM_1D_TE_RT_multi(epsilon, periodx, periody, h, lambda, theta, phi, refIndices, N, M, L)
%
% INPUT:
% <epsilon>: vector with permitivity distribution
% <period>: period of the grating in micron
% <lambda>: wavelength in micron
% <theta, phi>: angle of incidence in degree
% <refIndices>: vector [n1 n2] with refr.indices of surrounding
% <N>: number of Fourier orders


% OUTPUT:
% <gamma>: vector with propagation constants of grating modes
% program for TE polarization
% for TM polarization we switch eps <-> -mu and E<->H

%after Li,1997 + Li,Chapter 13 + Lecture 17,18(CEM)

k0 = 2*pi/lambda;
NN = (2*N+1)*(2*N+1);

gamma = zeros(4*NN,4*NN,L);

alpha_v = zeros(2*N+1,1);
beta_v = zeros(2*N+1,1);
alpha0 = refIndices(1)*k0*sin(theta)*cos(phi);
beta0 = refIndices(1)*k0*sin(theta)*sin(phi);
%gamma0 = -refIndices(1)*k0*cos(theta);  %or the same with "+"? Chapter 13, p.35.

for m=1:(2*N+1)
    alpha_v(m) = (alpha0 + (m-N-1)*2*pi/periodx);
    beta_v(m) = (beta0 + (m-N-1)*2*pi/periody);
end
alpha = zeros(NN,NN);
beta = zeros(NN,NN);
K = 2*N+1;
j=1;
%{
for i = 1:(2*N+1)*(2*N+1)
    if (floor((i-2)/K) < floor((i-1)/K)) && (i>2)
    j = j+1;
    end
    alpha(i,i) = alpha_v(j);
end
%}
for i = 1:(2*N+1)
    for j=1:(2*N+1)
        m = i+(2*N+1)*(j-1);
        beta(m,m) = beta_v(i);
        alpha(m,m) = alpha_v(j);       
    end
end

W = zeros(4*NN,4*NN,L);
pplus = zeros(2*NN,2*NN,L);
pminus = zeros(2*NN,2*NN,L);


for i=1:L
    [gammat, Wt, pplust, pminust, eps] = FMM_1D_TE_RT_beta_e(alpha, beta, epsilon(:,:,i), periodx, periody, h(i), lambda, theta, phi, refIndices, N, M);
    gamma(:,:,i) = gammat;
    W(:,:,i) = Wt;
    pplus(:,:,i) = pplust;
    pminus(:,:,i) = pminust;
end


kz1v = zeros(NN,1);
kz2v = zeros(NN,1);
A1 = zeros(NN,NN);
B1 = zeros(NN,NN);
C1 = zeros(NN,NN);
A2 = zeros(NN,NN);
B2 = zeros(NN,NN);
C2 = zeros(NN,NN);
k1 = refIndices(1)*k0;
k2 = refIndices(2)*k0;

for i = 1:(2*N+1)
    for j=1:(2*N+1)
        m = i+(2*N+1)*(j-1);
        kz1v(m) = ( k1^2 - (alpha_v(j))^2 - (beta_v(i))^2 )^(1/2);
        kz2v(m) = ( k2^2 - (alpha_v(j))^2 - (beta_v(i))^2 )^(1/2);
        A1(m,m) = ( k1^2 - (alpha_v(j))^2)/(k0*kz1v(m));
        A2(m,m) = ( k2^2 - (alpha_v(j))^2)/(k0*kz2v(m));
        B1(m,m) = ( k1^2 - (beta_v(i))^2)/(k0*kz1v(m));
        B2(m,m) = ( k2^2 - (beta_v(i))^2)/(k0*kz2v(m));
        C1(m,m) = alpha_v(j)*beta_v(i)/(k0*kz1v(m));
        C2(m,m) = alpha_v(j)*beta_v(i)/(k0*kz2v(m));
    end
end

kz1norm = kz1v/k0
mzero = zeros(NN,NN);
miden = eye(NN,NN);

K1_1 = cat(2, miden, mzero, miden, mzero);
K1_2 = cat(2, mzero, miden, mzero, miden);
K1_3 = cat(2, -C1, -A1, C1, A1);
K1_4 = cat(2, B1, C1, -B1, -C1);
K1 = cat(1, K1_1, K1_2, K1_3, K1_4);

K2_1 = cat(2, miden, mzero, miden, mzero);
K2_2 = cat(2, mzero, miden, mzero, miden);
K2_3 = cat(2, -C2, -A2, C2, A2);
K2_4 = cat(2, B2, C2, -B2, -C2);
K2 = cat(1, K2_1, K2_2, K2_3, K2_4);


%{
Smin1 = eye(4*NN,4*NN);
S0 = new_recursion(Smin1, K2, W(:,:,L), eye(2*NN,2*NN), eye(2*NN,2*NN), N);
Stemp = S0;
if L>1
    for i=(L-1):1
        Si = new_recursion(Stemp, W(:,:,i+1), W(:,:,i), pplus(:,:,i), pminus(:,:,i), N);
        Stemp = Si;
    end
end
Stotal = new_recursion(Stemp, W(:,:,1), K1, pplus(:,:,L), pminus(:,:,L), N);
%}


Smin1 = eye(4*NN,4*NN);
S0 = new_recursion_book(Smin1, K2, W(:,:,1), eye(2*NN,2*NN), pminus(:,:,1), N);
Stemp = S0;
if L>1
    for i=1:(L-1)
        Si = new_recursion_book(Stemp, W(:,:,i), W(:,:,i+1), pplus(:,:,i), pminus(:,:,i+1), N);
        Stemp = Si;
    end
end
Stotal = new_recursion_book(Stemp, W(:,:,L), K1, pplus(:,:,L), eye(2*NN,2*NN), N);

i = N+1;
j = N+1;
nul = i+(2*N+1)*(j-1);
%norm = A1(nul,nul) + B1(nul,nul) + C1(nul,nul)*2;
norm = A1(nul,nul);
I1 = 1/sqrt(norm);
I2 = 1/sqrt(norm);
u01 = zeros(NN,1);
u02 = zeros(NN,1);
       
dlast1 = zeros(NN,1);
dlast2 = zeros(NN,1);
dlast2(nul) = I2;
plast = zeros(2*NN,2*NN);
for i=1:NN
    %plast(2*i-1,2*i-1)=exp(-1i*kz1v(i)*h);
    %plast(2*i,2*i) = exp(-1i*kz1v(i)*h);
    plast(i,i) = exp(-1i*kz1v(i)*h);
    plast(i+NN,i+NN) = exp(-1i*kz1v(i)*h);
end
d = plast*cat(1,dlast1,dlast2);
delta = cat(1,u01,u02,d);
%delta = cat(1,dlast1,dlast2,u01,u02);

R_T = Stotal*delta;
%R = R_T(1:2*NN);
%T = R_T((2*NN+1):4*NN);
%T = pminus(:,:,L)\T;
%R_T=cat(1,R,T);

R1 = R_T(1:NN);        %u_last
R2 = R_T((NN+1):2*NN);
T1 = R_T((2*NN+1):3*NN);  %d0
T2 = R_T((3*NN+1):4*NN);

%nR = size(R1)
%nnT_R = size(T_R)
eta_R = zeros(NN,1);
eta_T = zeros(NN,1);

for i=1:NN
    eta_R(i) = A1(i,i)*(abs(R2(i)))^2 + B1(i,i)*(abs(R1(i)))^2 + C1(i,i)*( R1(i)*conj(R2(i))+R2(i)*conj(R1(i)) );
    eta_T(i) = A2(i,i)*(abs(T2(i)))^2 + B2(i,i)*(abs(T1(i)))^2 + C2(i,i)*( T1(i)*conj(T2(i))+T2(i)*conj(T1(i)) );
end


end