clc
clear all
load MatParam_nk_Au_2_interpExportData.txt
load MatParam_nk_ZnSe_interpExportData.txt
load MatParam_nk_SiO2_interpExportData.txt
load MatParam_nk_Zep520A_interpExportData.txt
profile on -history

lamnkAu = MatParam_nk_Au_2_interpExportData;
lamnkZnSe = MatParam_nk_ZnSe_interpExportData;
lamnkSiO2 = MatParam_nk_SiO2_interpExportData;
lamnkZep = MatParam_nk_Zep520A_interpExportData;

lmin = 600;
lmax = 1600;
tmin = 0;
tmax = 60;
lambda = linspace(lmin,lmax,45)*10^(-9);
lambdanm = linspace(lmin,lmax,45);
[Nll, Nl] = size(lambda)

%for i=1:Nl
%num = lambdanm(i) - 500 + 1;
num = floor(lambdanm - 500 + 1)
nAu = lamnkAu(num,2) + 1j*lamnkAu(num,3);
nZnSe = lamnkZnSe(num,2) + 1j*lamnkZnSe(num,3);
nSiO2 = lamnkSiO2(num,2) + 1j*lamnkSiO2(num,3);
nZep = lamnkZep(num,2) + 1j*lamnkZep(num,3);
%end

N = 4;                %number of Fourier orders
L = 3;                 %number of layers
periodx = 200*(10^(-9));  %period of periodic layer
periody = 200*(10^(-9));  %period of periodic layer
rx = (140/2)*(10^(-9));        %disc radius
ry = (175/2)*(10^(-9));        %disc radius
a = periodx;  
h = zeros(L,1);
h(3) = 25*10^(-9);
h(2) = 100*10^(-9);
h(1) = 20*10^(-9);

M = 301;               %number of modes for Fourier transform of epsilon
Mrx = (rx/a)*M;
Mry = (ry/a)*M;
K=1;
epsilon = zeros(M, M, L, Nl); %K - number of options
i0 = 1+floor(M/2);
j0 = 1+floor(M/2);
for k=1:Nl
    epsilon(:,:,3,k) = ones(M,M);
    epsilon(:,:,2,k) = nSiO2(k)*ones(M,M);
    epsilon(:,:,1,k) = ones(M,M);
    for i=1:(i0-1)
    for j=1:(j0-1)         
        if ( ((i-i0)^2/(Mrx^2)+(j-j0)^2/(Mry^2)) <= 1)    
        epsilon(j,i,1,k) = nAu(k)^2;
        epsilon(M-j+1,i,1,k) = nAu(k)^2;
        epsilon(j,M-i+1,1,k) = nAu(k)^2;
        epsilon(M-j+1,M-i+1,1,k) = nAu(k)^2;
        end
    end
    end
    i = i0;
    for j=1:M
        if ( (j-j0)^2/(Mry^2) <= 1)    
        epsilon(j,i,2,k) = nAu(k)^2;
        end
    end
    j = j0;
    for i=1:M
        if ( (i-i0)^2/(Mrx^2) <= 1)    
        epsilon(j,i,2,k) = nAu(k)^2;
        end
    end
end
refIndices = zeros(k,2);
for k=1:Nl
refIndices(k,:) = [nZnSe(k) 1.0];
end

%phi = linspace(0,90,40)*pi/180;
%[Npp,Np] = size(phi)
phi = 0*pi/180;
Np = 1;
theta = linspace(tmin,tmax,35)*pi/180;
[Ntt, Nt] = size(theta)

%Rsum = zeros(Nt,Np);
%Tsum = zeros(Nt,Np);
Rsum = zeros(Nl,Nt);
Tsum = zeros(Nl,Nt);

P = 2*N+1;
Q = 2*N+1;
R = 1;

eps11f=zeros(P*Q,P*Q,L,Nl);
eps22f=zeros(P*Q,P*Q,L,Nl);
eps33f=zeros(P*Q,P*Q,L,Nl);
[eps11air(:,:), eps22air(:,:), eps33air(:,:)] = FMM_eps123_new(epsilon(:,:,3,1),N,M);


for k=1:Nl
    for i=1:(L-1)
    [eps11f(:,:,i,k), eps22f(:,:,i,k), eps33f(:,:,i,k)] = FMM_eps123_new(epsilon(:,:,i,k),N,M);
    end
    eps11f(:,:,3,k) = eps11air;
    eps22f(:,:,3,k) = eps22air;
    eps33f(:,:,3,k) = eps33air;
        
end



for i=1:Nl
    for j=1:Nt
    for k=1:Np 
    eps11 = eps11f(:,:,:,i);
    eps22 = eps22f(:,:,:,i);
    eps33 = eps33f(:,:,:,i);
    [eta_R, eta_T] = FMM_1D_TE_RT_multi(eps11,eps22,eps33, periodx,periody, h, lambda(i), theta(j), phi(k), refIndices(i,:), N, M, L);
    %[gamma, eta_R, eta_T] = FMM_1D_TE_RT_multi(eps11d,eps22d,eps33d,epsilon, periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L);
    Rsum(i,j) = sum(eta_R);
    Tsum(i,j) = sum(eta_T);
    %Rsum(j,k) = sum(eta_R);
    %Tsum(j,k) = sum(eta_T);
    end
    end
    ii=i
end
%{
k0=2*pi/lambda;
ka = 2*nZnSe*a/lambda;
kx = ka*sin(theta')*cos(phi);
ky = ka*sin(theta')*sin(phi);
figure;
pcolor(kx,ky,Rsum)
shading flat
caxis([0 1])
%}
[l, ll] = size(lambda)
[t,tt] = size(theta)
[r,rr] = size(Rsum)
tl = linspace(lmin,lmax,200)/1000;
tt = linspace(tmin,tmax,200);
[XI,YI] = meshgrid(tl,tt);
ZI = griddata(lambdanm/1000,theta*180/pi,transpose(Rsum),XI,YI);
figure;
pcolor(XI,YI,ZI)
shading flat
caxis([0 1])
%{
figure(2)
hold on
plot(lambda, Tsum./Tsum0)
hold off
%}
p = profile('info')

