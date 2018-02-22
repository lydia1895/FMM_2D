clc
clear all

N = 3;                %number of Fourier orders
L = 3; %number of layers
%periodx = zeros(L,1);
%periody = zeros(L,1);
wx = zeros(L,1);
wy = zeros(L,1);
periodx_top_real = 1300*10^(-9);
periodx = 3900*10^(-9);
periody = 3000*10^(-9);
wx(3) = 1000*10^(-9);
wy(2) = 1000*10^(-9);
wx(1) = 1000*10^(-9);

h = zeros(L,1);
h(3) = 1000*10^(-9);
h(2) = 1000*10^(-9);
h(1) = 1000*10^(-9);
M = 10001;

epsilon = zeros(M, M, L);


nSiO2 = 1.45;
epsSiO2 = nSiO2^2;

n_resist = 1.53;
eps_resist = n_resist^2;

n_air = 1.0;
eps_air = 1.0;

epsilon(:,:,3)=eps_air*ones(M,M,1);
epsilon(:,:,2)=eps_air*ones(M,M,1);
epsilon(:,:,1)=eps_air*ones(M,M,1);
refIndices = [n_air n_air];

%in epsilon row is y, column is x
%lowest layer has nlayer=1, top layer has nlayer=L
x = zeros(M,L);
y = zeros(M,L);
x(:,3) = (1:1:M)*periodx/M;
y(:,2) = (1:1:M)*periody/M;
x(:,1) = (1:1:M)*periodx/M;
for i=1:M
    for j=1:M
        if (x(i,3)<wx(3)) || (periodx_top_real<x(i,3))&&(x(i,3)<periodx_top_real+wx(3)) ||...
                (2*periodx_top_real<x(i,3))&&(x(i,3)<2*periodx_top_real+wx(3))
            epsilon(j,i,3) = eps_resist;
        end
        if y(j,2)<wy(2)
            epsilon(j,i,2) = eps_resist;
        end
        if x(i,1)<wx(1)
            epsilon(j,i,1) = eps_resist;
        end
    end
end

lmin = 1400*10^(-9);
lmax = 1550*10^(-9);
thetamin = 0;
thetamax = 1*pi/180;
lambda = linspace(lmin, lmax, 200);
[Nll,Nl] = size(lambda);
%theta = linspace(thetamin,thetamax,20);
%[Ntt,Nt]=size(theta);
theta = [0.1 0.5 1]*pi/180;
Nt=3;

phi = 0*pi/180;
Np=1;
Rsum=zeros(Nl,Nt);
Tsum=zeros(Nl,Nt);

polarization = 'TE';

P = 2*N+1;
Q = 2*N+1;
R = 1;

eps11=zeros(P*Q,P*Q,L);
eps22=zeros(P*Q,P*Q,L);
eps33=zeros(P*Q,P*Q,L);
for i=1:L
    [eps11(:,:,i), eps22(:,:,i), eps33(:,:,i)] = FMM_eps123_new(epsilon(:,:,i),N,M);
end



for i=1:Nl
    for j=1:Nt
        for k=1:Np
            [eta_R1, eta_T1] = FMM_2D_RT_multi(eps11,eps22,eps33,...
                periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L,polarization);
            %Rsum(j,k) = sum(eta_R1);
            %Tsum(j,k) = sum(eta_T1);
            Rsum(i,j) = sum(eta_R1);
            Tsum(i,j) = sum(eta_T1);
        end
    end
end
%{
epsilon(:,:,2)=epsTa2O5*ones(M,M,1);

Rsum_non=zeros(Nl,Nt);
Tsum_non=zeros(Nl,Nt);

eps11=zeros(P*Q,P*Q,L);
eps22=zeros(P*Q,P*Q,L);
eps33=zeros(P*Q,P*Q,L);
for i=1:L
[eps11(:,:,i), eps22(:,:,i), eps33(:,:,i)] = FMM_eps123_new(epsilon(:,:,i),N,M);
end



for i=1:Nl
    for j=1:Nt
    for k=1:Np
    [eta_R1, eta_T1] = FMM_2D_RT_multi(eps11,eps22,eps33,...
        periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L,polarization);
    %Rsum(j,k) = sum(eta_R1);
    %Tsum(j,k) = sum(eta_T1);
    Rsum_non(i,j) = sum(eta_R1);
    Tsum_non(i,j) = sum(eta_T1);
    end
    end
end
%}

lambda = lambda*10^6;
lmin = lmin*10^6;
lmax = lmax*10^6;



figure(1)
plot(lambda, Rsum(:,1), 'b', lambda, Rsum(:,2)+0.5, 'g', lambda, Rsum(:,3)+1.0,'r')%,...
    %lambda, Rsum(:,4)+1.5, 'm', 'LineWidth', 2)
h5 = legend('theta=0.1 deg','theta=0.5 deg','theta=1 deg',3);
set(h5,'Interpreter','none')
axis tight
ax = gca;
ax.XAxis.MinorTick = 'on';
xlabel('lambda, mkm')
ylabel('R')
set(gca,'fontsize', 16)



%{
figure(1)
tl = linspace(lmin,lmax,400)*10^9;
tt = linspace(thetamin, thetamax,40)*180/pi;
[XI,YI] = meshgrid(tl,tt);
ZI = griddata(theta*180/pi,lambda*10^9,Rsum,YI,XI);
figure;
%pcolor(kx,ky,Rsum)
pcolor(XI,YI,ZI)
shading flat
caxis([0 1])
colormap('jet');
colorbar;
hold off
%}
%{
k0=2*pi/lambda;
kx = sin(theta')*cos(phi);
ky = sin(theta')*sin(phi);
figure;
pcolor(kx,ky,Rsum)
%shading interp
colormap gray;
%}
%pcolor(theta*180/pi,lambda,Rsum);
%imagesc(theta*180/pi,lambda,Rsum);
%set(gca,'Yscale','linear','Ydir','normal');

