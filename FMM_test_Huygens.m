clc
clear all
profile on -history


N = 3;                %number of Fourier orders
L = 2;                 %number of layers
periodx = 666*(10^(-9));  %period of periodic layer
r = 242*(10^(-9));        %disc radius
periody = 666*(10^(-9));  %period of periodic layer
a = periodx;  
h = zeros(L,1);
h(2) = 220*(10^(-9));       %thickness of periodic layer
h(1) = 2*10^(-6);
%h(1) = 10^(-3);
M = 301;               %number of modes for Fourier transform of epsilon
Mr = (r/a)*M;
K = 2;
epsilon = zeros(M, M, L, K); %number of options
i0 = 1+floor(M/2);
j0 = 1+floor(M/2);
for i=1:M
    for j=1:M
    %epsilon(i,j,1,1) = 3.5^2;
    %epsilon(i,j,1,2) = 3.5^2;    
    epsilon(i,j,1,1) = 1.45^2;
    epsilon(i,j,1,2) = 1.45^2;
    if ( ((i-i0)^2+(j-j0)^2) <= Mr^2)    
        epsilon(j,i,2,1) = 1.0^2;
        epsilon(j,i,2,2) = 3.5^2;
    else
        epsilon(j,i,2,1) = 1.0^2;
        epsilon(j,i,2,2) = 1.0^2;
    end
    end
end
refIndices = [1.0 3.5];


lambda = linspace(1000,1700,50)*10^(-9);
[Nll,Nl] = size(lambda)
theta = 0*pi/180;
Nt=1;
phi = 0*pi/180;
Np=1;
Rsum=zeros(Nl,Nt);
Tsum=zeros(Nl,Nt);
Rsum0=zeros(Nl,Nt);
Tsum0=zeros(Nl,Nt);


P = 2*N+1;
Q = 2*N+1;
R = 1;

eps11=zeros(P*Q,P*Q,L,K);
eps22=zeros(P*Q,P*Q,L,K);
eps33=zeros(P*Q,P*Q,L,K);
for i=1:L
    for k=1:K
    [eps11(:,:,i,k), eps22(:,:,i,k), eps33(:,:,i,k)] = FMM_eps123_new(epsilon(:,:,i,k),N,M);
    end
end
eps110 = eps11(:,:,:,1);
eps220 = eps22(:,:,:,1);
eps330 = eps33(:,:,:,1);
eps11d = eps11(:,:,:,2);
eps22d = eps22(:,:,:,2);
eps33d = eps33(:,:,:,2);

for i=1:Nl
    for j=1:Nt
    for k=1:Np 
    [gamma0, eta_R0, eta_T0] = FMM_1D_TE_RT_multi(eps110,eps220,eps330,epsilon, periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L);
    [gamma, eta_R, eta_T] = FMM_1D_TE_RT_multi(eps11d,eps22d,eps33d,epsilon, periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L);
    Rsum0(i,j) = sum(eta_R0);
    Tsum0(i,j) = sum(eta_T0);
    Rsum(i,j) = sum(eta_R);
    Tsum(i,j) = sum(eta_T);
    end
    end
end

figure(2)
hold on
plot(lambda, Tsum./Tsum0)
hold off

p = profile('info')
