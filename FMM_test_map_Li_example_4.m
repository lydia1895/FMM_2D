clc
clear all
load MatParam_nk_Au_2_interpExportData.txt
load MatParam_nk_ZnSe_interpExportData.txt
load MatParam_nk_SiO2_interpExportData.txt
load MatParam_nk_Zep520A_interpExportData.txt
%profile on -history


lamnkAu = MatParam_nk_Au_2_interpExportData;
lamnkZnSe = MatParam_nk_ZnSe_interpExportData;
lamnkSiO2 = MatParam_nk_SiO2_interpExportData;
lamnkZep = MatParam_nk_Zep520A_interpExportData;

%lambda = 950*(10^(-9));
%lambdanm = 950;
lmin = 600;
lmax = 1300;
lambda = linspace(lmin,lmax,40)*10^(-9);
lambdanm = linspace(lmin,lmax,40);
[Nll, Nl] = size(lambda)

nAu = zeros(Nl,1);
nZnSe = zeros(Nl,1);
nSiO2 = zeros(Nl,1);
nZep = zeros(Nl,1);
for i=1:Nl
num = floor(lambdanm(i)) - 500 + 1;
nAu(i) = lamnkAu(num,2) + 1j*lamnkAu(num,3)
nZnSe(i) = lamnkZnSe(num,2) + 1j*lamnkZnSe(num,3);
nSiO2(i) = lamnkSiO2(num,2) + 1j*lamnkSiO2(num,3);
nZep(i) = lamnkZep(num,2) + 1j*lamnkZep(num,3);
end

N = 5;                %number of Fourier orders
L = 3;                 %number of layers
periodx = 200*(10^(-9));  %period of periodic layer
periody = 200*(10^(-9));  %period of periodic layer
rx = (140/2)*(10^(-9));        %disc radius
ry = (175/2)*(10^(-9));        %disc radius
a = periodx;  
h = zeros(L,1);
h(3) = 25*10^(-9);
h(2) = 180*10^(-9);       %thickness of periodic layer
h(1) = 20*10^(-9);
M = 201;               %number of modes for Fourier transform of epsilon
Mrx = (rx/a)*M;
Mry = (ry/a)*M;

epsilon = zeros(M, M, L, Nl); %K - number of options
refIndices = zeros(Nl,2);
i0 = 1+floor(M/2);
j0 = 1+floor(M/2);
for k =1:Nl    
    for i=1:M
    for j=1:M   
    epsilon(j,i,3,k) = 1.0;
    epsilon(j,i,2,k) = nZep(k)^2;
    if ( ((i-i0)^2/(Mrx^2)+(j-j0)^2/(Mry^2)) <= 1)    
        epsilon(j,i,1,k) = nAu(k)^2;
    else
        epsilon(j,i,1,k) = nZep(k)^2;
    end
    end
    end
refIndices(k,:) = [nZnSe(k) nSiO2(k)];
end




%phi = linspace(0,90,40)*pi/180;
%[Npp,Np] = size(phi)
%tmin = 0;
%tmax = 60;
%theta = linspace(tmin,tmax,30)*pi/180;
%[Ntt, Nt] = size(theta);
theta = 40*pi/180;
Nt=1;
%theta = [40 55]*pi/180;
%Nt=2;
phi = 45*pi/180;
Np = 1;


Rsum = zeros(Nl,Nt);
Tsum = zeros(Nl,Nt);
%Nl=1;
%Rsum = zeros(Nt,Np);
%Tsum = zeros(Nt,Np);

P = 2*N+1;
Q = 2*N+1;
R = 1;

eps11f=zeros(P*Q,P*Q,L,Nl);
eps22f=zeros(P*Q,P*Q,L,Nl);
eps33f=zeros(P*Q,P*Q,L,Nl);
for i=1:L
    for k=1:Nl
    [eps11f(:,:,i,k), eps22f(:,:,i,k), eps33f(:,:,i,k)] = FMM_eps123_new(epsilon(:,:,i,k),N,M);
    epsk = k
    end
end



for i=1:Nl
    for j=1:Nt
    for k=1:Np
    eps11 = eps33f(:,:,:,i);
    eps22 = eps33f(:,:,:,i);
    eps33 = eps33f(:,:,:,i); 
    [eta_R, eta_T] = FMM_1D_TE_RT_multi(eps11,eps22,eps33,periodx, periody, h, lambda(i), theta(j), phi(k), refIndices(i,:), N, M, L);
    %[eta_R, eta_T] = FMM_1D_TE_RT_multi(eps11d,eps22d,eps33d,epsilon, periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L);
    Rsum(i,j) = sum(eta_R);
    Tsum(i,j) = sum(eta_T);
    ii = i
    jj = j
    %Rsum(j,k) = sum(eta_R);
    %Tsum(j,k) = sum(eta_T);
    
    end
    end
end
%{
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
%}
%{
k0=2*pi/lambda;
ka = 2*nZnSe*a/lambda;
kx = ka*sin(theta')*cos(phi);
ky = ka*sin(theta')*sin(phi);

ti = linspace(0,ka,100);
[XI,YI] = meshgrid(ti,ti);
ZI = griddata(kx,ky,Rsum,XI,YI);
figure;
%pcolor(kx,ky,Rsum)
pcolor(XI,YI,ZI)
shading flat
caxis([0 1])
%}


figure(2)
hold on
plot(lambda, Rsum(:,1), 'b')%, lambda, Rsum(:,2), 'r')
hold off

%p = profile('info')

