function [eta_R, eta_T] = FMM_1D_TE_RT_multi(eps11,eps22,eps33,periodx, periody, h, lambda, theta, phi, refIndices, N, M, L)
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
%gamma01 = -refIndices(1)*k0*cos(theta);  %or the same with "+"? Chapter 13, p.35.
k1 = refIndices(1)*k0;
k2 = refIndices(2)*k0;
gamma01 = ( k1^2 - (alpha0)^2 - (beta0)^2 )^(1/2);
gamma02 = ( k2^2 - (alpha0)^2 - (beta0)^2 )^(1/2);

for m=1:(2*N+1)
    alpha_v(m) = (alpha0 + (m-N-1)*2*pi/periodx);
    beta_v(m) = (beta0 + (m-N-1)*2*pi/periody);
end
alpha = zeros(NN,NN);
beta = zeros(NN,NN);

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
    [Wt, pplust, pminust] = FMM_2D_gamma_e(eps11(:,:,i),eps22(:,:,i),eps33(:,:,i),alpha, beta, h(i),lambda,N);
    %[gammat, Wt, pplust, pminust, eps] = FMM_1D_TE_RT_beta_e(alpha, beta,
    %epsilon(:,:,i), periodx, periody, h(i), lambda, theta, phi, refIndices, N, M);
    W(:,:,i) = Wt;
    pplus(:,:,i) = pplust;
    pminus(:,:,i) = pminust;
end
%eps22=0;
%eps33=0;

kz1v = zeros(NN,1);
kz2v = zeros(NN,1);
A1 = zeros(NN,NN);
B1 = zeros(NN,NN);
C1 = zeros(NN,NN);
A2 = zeros(NN,NN);
B2 = zeros(NN,NN);
C2 = zeros(NN,NN);


for i = 1:(2*N+1)
    for j=1:(2*N+1)
        m = i+(2*N+1)*(j-1);
        %beta(i+(2*N+1)*(j-1),i+(2*N+1)*(j-1)) = beta_v(i);    
        
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

kz1norm = kz1v/k0;
mzero = zeros(NN,NN);
miden = eye(NN,NN);

K1_1 = cat(2, miden, mzero, miden, mzero);
K1_2 = cat(2, mzero, miden, mzero, miden);
K1_3 = cat(2, -C1, -A1, C1, A1);
%K1_3 = cat(2, -C1, A1, C1, -A1);
K1_4 = cat(2, B1, C1, -B1, -C1);
K1 = cat(1, K1_1, K1_2, K1_3, K1_4);   %'2' in Chapter 13, Li

K2_1 = cat(2, miden, mzero, miden, mzero);
K2_2 = cat(2, mzero, miden, mzero, miden);
K2_3 = cat(2, -C2, -A2, C2, A2);
%K2_3 = cat(2, -C2, A2, C2, -A2);
K2_4 = cat(2, B2, C2, -B2, -C2);
K2 = cat(1, K2_1, K2_2, K2_3, K2_4);   %'0' in Chapter 13, Li



%{
Smin1 = eye(4*NN,4*NN);
S0 = new_recursion(Smin1, K2, W(:,:,1), eye(2*NN,2*NN), pminus(:,:,1), N);
Stemp = S0;
if L>1
    for i=1:(L-1)
        Si = new_recursion(Stemp, W(:,:,i), W(:,:,i+1), pplus(:,:,i), pminus(:,:,i+1), N);
        Stemp = Si;
    end
end
Stotal = new_recursion(Stemp, W(:,:,L), K1, pplus(:,:,L), eye(2*NN,2*NN), N);
%}

Smin1 = eye(4*NN,4*NN);
Rudmin1 = zeros(2*NN,2*NN);
Rud0 = new_recursion_refl_only(Rudmin1, K2, W(:,:,1), eye(2*NN,2*NN), pminus(:,:,1), N);
Rudtemp = Rud0;
if L>1
    for i=1:(L-1)
        Rudi = new_recursion_refl_only(Rudtemp, W(:,:,i), W(:,:,i+1), pplus(:,:,i), pminus(:,:,i+1), N);
        Rudtemp = Rudi;
    end
end
Rudtotal = new_recursion_refl_only(Rudtemp, W(:,:,L), K1, pplus(:,:,L), eye(2*NN,2*NN), N);


i = N+1;
j = N+1;
nul = i+(2*N+1)*(j-1);
%TE modes
I1 = -sin(phi);
I2 = cos(phi);
%TM modes
%I1 = cos(phi)*cos(theta);
%I2 = sin(phi)*cos(theta);
norm = A1(nul,nul)*abs(I2)^2 + B1(nul,nul)*abs(I1)^2 + C1(nul,nul)*( I1*conj(I2)+I2*conj(I1) );
I1=I1/sqrt(norm);
I2=I2/sqrt(norm);

u01 = zeros(NN,1);
u02 = zeros(NN,1);
       
dlast1 = zeros(NN,1);
dlast2 = zeros(NN,1);
dlast1(nul) = I1*exp(1j*gamma01*sum(h));
dlast2(nul) = I2*exp(1j*gamma01*sum(h));
delta = cat(1,u01,u02,dlast1,dlast2);

dlast = cat(1,dlast1,dlast2);
u = Rudtotal*dlast;
R1 = zeros(NN);
R2 = zeros(NN);
for i=1:NN
    R1(i) = u(i)/exp(1j*kz1v(i)*(sum(h)-h(L)));
    R2(i) = u(i+NN)/exp(1j*kz1v(i)*(sum(h)-h(L)));
end
%R1 = R(1:NN)/exp( 1j*gamma01*(sum(h)-h(L)) );        %u_last
%R2 = R((NN+1):2*NN)/exp( 1j*gamma01*(sum(h)-h(L)) );

%{
R_T = Stotal*delta;
R1 = R_T(1:NN)/exp( 1j*gamma01*(sum(h)-h(L)) );        %u_last
R2 = R_T((NN+1):2*NN)/exp( 1j*gamma01*(sum(h)-h(L)) );
T1 = R_T((2*NN+1):3*NN);  %d0
T2 = R_T((3*NN+1):4*NN);
%}

eta_R = zeros(NN,1);
eta_T = zeros(NN,1);

for i=1:NN
    if imag(kz1v(i))==0
    eta_R(i) = A1(i,i)*(abs(R2(i)))^2 + B1(i,i)*(abs(R1(i)))^2 + C1(i,i)*( R1(i)*conj(R2(i))+R2(i)*conj(R1(i)) );
    end
    %if imag(kz2v(i))==0
    %eta_T(i) = A2(i,i)*(abs(T2(i)))^2 + B2(i,i)*(abs(T1(i)))^2 + C2(i,i)*( T1(i)*conj(T2(i))+T2(i)*conj(T1(i)) );
    %end
end

end