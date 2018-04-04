clc
clear all


N = 5;                %number of Fourier orders
L = 1;                 %number of layers

%?????? 450, 500, 550 ??. ?????? ??????? 500 ??, ?????? 300 ??.
%?????? ???? ????? 100 ??. 
%? ???????? ???????? ????? ?????? ITO ? ???? ?????????? ??????????? ???????? 1.9.
%???????? ???? ???? ????????? 530-560 ??.

periodx = 550;  %period of periodic layer
periody = 3000;  %period of periodic layer
wx_etched = 100;
wx_width = periodx - wx_etched;
wy_width = 500;
h = zeros(L,1);
h(2) = 250;
h(1) = 300;       %thickness of periodic layer

M = 90;               %number of modes for Fourier transform of epsilon
x = (1:1:M)*periodx/M;
y = (1:1:M)*periody/M;

lmin = 501;
lmax = 700;
Nl=200;
lambda = linspace(lmin,lmax,Nl);

n_air = 1.0;
eps_air = n_air^2;
%n_prism = 2.3*ones(Nl,1);

n_Perovskite = zeros(Nl,1);
eps_Perovskite = zeros(Nl,1);
n_ITO = zeros(Nl,1);
eps_ITO = zeros(Nl,1);

load S1.txt
load ITO.txt

Perovskite_lambda = S1(:,1);
ITO_lambda = ITO(:,1)*1000;

for i=1:Nl
    [ll,num1] = min( abs (lambda(i)-Perovskite_lambda(:) ) );
    n_Perovskite(i) = S1(num1,2)+ 1j*S1(num1,3);
    eps_Perovskite(i) = n_Perovskite(i)^2;
    
    [ll,num2] = min( abs (lambda(i)-ITO_lambda(:) ) );
    n_ITO(i) = ITO(num2,2);%+ 1j*ITO(num2,3);
    eps_ITO(i) = n_ITO(i)^2;
end

%{
thetamin = 0*pi/180;
thetamax = 80*pi/180;
Nt=41;
theta = linspace(thetamin,thetamax,Nt);
%}

theta = 0*pi/180;
Nt=1;

%{
theta = [0.1 1 3]*pi/180;
Nt = 3;
%}
phi = 0*pi/180;
Np=1;
Rsum=zeros(Nl,Nt);
Tsum=zeros(Nl,Nt);


P = 2*N+1;
Q = 2*N+1;
R = 1;

eps11=zeros(P*Q,P*Q,L);
eps22=zeros(P*Q,P*Q,L);
eps33=zeros(P*Q,P*Q,L);


epsilon = zeros(M, M, L);
nn_air = ones(Nl,1);
refIndices = [nn_air n_ITO];
n_substrate = n_ITO;


%phase_R = zeros(Nl,Nt);
%Rsum_full = zeros(Nl,Nt);
polarization = 'TE';
%{
nlayer = 2;
epsilon(:,:,nlayer) = eps_air*ones(M,M);
[eps11(:,:,nlayer), eps22(:,:,nlayer), eps33(:,:,nlayer)] =...
        FMM_eps123_new(epsilon(:,:,nlayer),N,M);
%}
for i=1:Nl
    
    for ii=1:M
        for jj=1:M
            if ( x(ii)<=wx_width ) && ( y(jj)<=wy_width )
                epsilon(jj,ii,1) = eps_Perovskite(i);
            else
                epsilon(jj,ii,1) = eps_air;
                
            end
        end
    end
    nlayer = 1;
    [eps11(:,:,nlayer), eps22(:,:,nlayer), eps33(:,:,nlayer)] =...
        FMM_eps123_new(epsilon(:,:,nlayer),N,M);
    rrefIndices = [refIndices(i,1) refIndices(i,2)]  
    for j=1:Nt
        for k=1:Np
            
            [eta_R, eta_T, eta_R_full, eta_T_full] = FMM_1D_TE_RT_multi(eps11,eps22,eps33,...
                periodx, periody, h, lambda(i), theta(j), phi(k), rrefIndices, N, M, L, polarization);
            
            Rsum(i,j) = sum(eta_R);
            Tsum(i,j) = sum(eta_T);
            Rsum_full(i,j) = sum(eta_R_full);
            phase_R(i,j) = angle(Rsum_full(i,j));
        end
    end
    lambda(i)
      
end
%{
figure(1)
plot(lambda, Rsum(:,1), 'b', lambda, Rsum(:,2)+0.01, 'g', lambda, Rsum(:,3)+0.02,'r',...
    'LineWidth', 2)
h5 = legend('theta=0.1 deg','theta=1 deg','theta=3 deg', 3);
set(h5,'Interpreter','none')
axis tight
ax = gca;
ax.XAxis.MinorTick = 'on';
xlabel('lambda, nm')
ylabel('R')
set(gca,'fontsize', 16)

hold off
%}

figure(2)
plot(lambda, 1-Rsum-Tsum, 'r', 'LineWidth', 2)
axis tight
ax = gca;
ax.XAxis.MinorTick = 'on';
xlabel('lambda, nm')
ylabel('absorption')
set(gca,'fontsize', 16)

hold off

%{
figure(1);
pcolor(lambda/1000,theta*180/pi,transpose(Rsum))

xlabel('lambda, mkm');
ylabel('theta, deg');
colormap('jet');
colorbar;
set(gca,'fontsize', 16)
shading flat
caxis([0 1])
colorbar
hold off

a = periodx;

c = physconst('LightSpeed');
h = 4.135666 * 10^(-15);
kx = (2*a*n_prism(1)./lambda)'*sin(theta);
frequency = (c*10^(-3)./lambda');

figure(2);
pcolor(kx,frequency,Rsum)
title('reflection')
xlabel('kx, \pi/a');
ylabel('frequency, THz');
colormap('jet');
colorbar;
set(gca,'fontsize', 16)
shading flat
caxis([0 1])
colorbar
hold on



plot(frequency*n_prism(1)*a*2*10^3/c,frequency,'b',...
    frequency*n_substrate(Nl)*a*2*10^3/c,frequency,'k','Linewidth',4)
hold off

figure(3);
pcolor(kx,frequency,phase_R*180/pi)
title('phase_R')
xlabel('kx, \pi/a');
ylabel('frequency, THz');
colormap('jet');
colorbar;
set(gca,'fontsize', 16)
shading flat
%caxis([-0.4 0])
hcb=colorbar
title(hcb,'phase, deg')
hold on

plot(frequency*n_prism(1)*a*2*10^3/c,frequency,'b',...
    frequency*n_substrate(1)*a*2*10^3/c,frequency,'k','Linewidth',4)
hold off
%}
