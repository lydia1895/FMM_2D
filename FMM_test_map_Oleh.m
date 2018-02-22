clc
clear all
load MatParam_nk_Au_2_interpExportData.txt
load MatParam_nk_ZnSe_interpExportData.txt
load MatParam_nk_SiO2_interpExportData.txt
load MatParam_nk_Zep520A_interpExportData.txt
load GOLD_McPeak_3001500_nm.txt
profile on -history

%{
lamnkAu = MatParam_nk_Au_2_interpExportData;
lamnkZnSe = MatParam_nk_ZnSe_interpExportData;
lamnkSiO2 = MatParam_nk_SiO2_interpExportData;
lamnkZep = MatParam_nk_Zep520A_interpExportData;
%}
frepsAu = GOLD_McPeak_3001500_nm;
[s1Au, s2Au] = size(frepsAu);

lambda = 730*(10^(-9));
lambdanm = 730;
lamAu = zeros(s1Au,1);
const = 299*10^3;
del0 = abs((const/frepsAu(1,1)) - lambdanm);
k0 = 1;
for k=1:s1Au
    lamAu(k) = (const/frepsAu(k,1));
    del = abs(lamAu(k)-lambdanm);
    if del<del0
        del0 = del;
        k0 = k;
    end
end  
epsAu = frepsAu(k0,2)+1j*frepsAu(k0,3)
%{
num = lambdanm - 500 + 1
nAu = lamnkAu(num,2) + 1j*lamnkAu(num,3);
nZnSe = lamnkZnSe(num,2) + 1j*lamnkZnSe(num,3);
nSiO2 = lamnkSiO2(num,2) + 1j*lamnkSiO2(num,3);
nZep = lamnkZep(num,2) + 1j*lamnkZep(num,3);
%}

N = 2;                %number of Fourier orders
L = 2;                 %number of layers
periodx = 200*(10^(-9));  %period of periodic layer
periody = 200*(10^(-9));  %period of periodic layer
rx = (134/2)*(10^(-9));        %disc radius
ry = (103/2)*(10^(-9));        %disc radius
a = periodx;  
h = zeros(L,1);
%h(3) = 25*10^(-9);
h(2) = 180*10^(-9);       %thickness of periodic layer
h(1) = 20*10^(-9);
M = 301;               %number of modes for Fourier transform of epsilon
Mrx = (rx/a)*M;
Mry = (ry/a)*M;
K=1;
epsilon = zeros(M, M, L, K); %K - number of options
i0 = 1+floor(M/2);
j0 = 1+floor(M/2);
epsilon(:,:,2,1) = 1.45^2*ones(M,M);
for i=1:M
    for j=1:M   
    %epsilon(j,i,3,1) = 1.0;
    %epsilon(j,i,2,1) = 2.1;
    if ( ((i-i0)^2/(Mrx^2)+(j-j0)^2/(Mry^2)) <= 1)    
        epsilon(j,i,1,1) = epsAu;
    else
        epsilon(j,i,1,1) = 1.45^2;
    end
    end
end
nZnSe = 2.48;
ndiel = 1.45;
refIndices = [nZnSe ndiel];


phi = linspace(0,90,35)*pi/180;
[Npp,Np] = size(phi);
theta0 = asin(ndiel/nZnSe)*180/pi;
theta = linspace(theta0,89,35)*pi/180;
[Ntt, Nt] = size(theta)

Nl = 1;
Rsum = zeros(Nt,Np);
Tsum = zeros(Nt,Np);

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
eps11 = eps11(:,:,:,1);
eps22 = eps22(:,:,:,1);
eps33 = eps33(:,:,:,1);


for i=1:Nl
    for j=1:Nt
    for k=1:Np 
    [eta_R, eta_T] = FMM_1D_TE_RT_multi(eps11,eps22,eps33,periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L);
    %[gamma, eta_R, eta_T] = FMM_1D_TE_RT_multi(eps11d,eps22d,eps33d,epsilon, periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L);
    %Rsum0(i,j) = sum(eta_R0);
    %Tsum0(i,j) = sum(eta_T0);
    Rsum(j,k) = sum(eta_R);
    Tsum(j,k) = sum(eta_T);
    jj = j
    end
    end
end
%{
thetanul = linspace(0, theta0-0.1, 35)*pi/180;
[Ntt0, Nt0] = size(thetanul); 
Rnul = zeros(Nt0, Np);
ttheta = cat(1,thetanul,theta);
RRsum = cat(1,Rnul,Rsum);
%}

kx = nZnSe*sin(theta')*cos(phi);
ky = nZnSe*sin(theta')*sin(phi);
%{
kx = nZnSe*sin(ttheta')*cos(phi);
ky = nZnSe*sin(ttheta')*sin(phi);
%}

figure;
pcolor(kx,ky,Rsum)
shading flat
caxis([0 1])

%{
ti = linspace(0,nZnSe,100);
[XI,YI] = meshgrid(ti,ti);
ZI = griddata(kx,ky,Rsum,XI,YI);
figure;
%pcolor(kx,ky,Rsum)
pcolor(XI,YI,ZI)
shading flat
caxis([0 1])
%}
%{
figure(2)
hold on
plot(lambda, Tsum./Tsum0)
hold off
%}
p = profile('info')

