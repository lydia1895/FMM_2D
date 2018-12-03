clc
clear all


N = 3;                %number of Fourier orders
L = 1;                 %number of layers
periodx = 666;  %period of periodic layer
r = 242;        %disc radius
periody = 666;  %period of periodic layer
a = periodx;  
h = zeros(L,1);
h(1) = 220;       %thickness of periodic layer

M = 301;               %number of modes for Fourier transform of epsilon
Mr = (r/a)*M;

epsilon = zeros(M, M, L); %number of options
i0 = 1+floor(M/2);
j0 = 1+floor(M/2);
for i=1:M
    for j=1:M
    if ( ((i-i0)^2+(j-j0)^2) <= Mr^2)    
        epsilon(j,i,1) = 3.5^2;      
    else
        epsilon(j,i,1) = 1.66^2;
        
    end
    end
end
refIndices = [1.66 1.66];

%lambda = 1350;
lambda = linspace(1300,1400,10);
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

eps11=zeros(P*Q,P*Q,L);
eps22=zeros(P*Q,P*Q,L);
eps33=zeros(P*Q,P*Q,L);
for i=1:L
    [eps11(:,:,i), eps22(:,:,i), eps33(:,:,i)] = FMM_eps123_new(epsilon(:,:,i),N,M);
end

polarization = 'TM';
for i=1:Nl
    for j=1:Nt
    for k=1:Np 
    %[gamma0, eta_R0, eta_T0] = FMM_1D_TE_RT_multi(eps110,eps220,eps330,epsilon, periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L);
    [eta_R, eta_T, eta_R_full, eta_T_full] = FMM_1D_TE_RT_multi(eps11,eps22,eps33,periodx, periody, h,...
        lambda(i), theta(j), phi(k), refIndices, N, M, L,polarization);
   
    Rsum(i,j) = sum(eta_R);
    Tsum(i,j) = sum(eta_T);
    Rsum_full(i,j) = sum(eta_R_full);
    Tsum_full(i,j) = sum(eta_T_full);
    end
    end
end

figure(1)
hold on
plot(lambda, Tsum,'r', lambda, angle(Tsum_full),'g')
hold off
figure(2)
hold on
plot(lambda, Rsum,'b', lambda, angle(Rsum_full),'g')
hold off
p = profile('info')
