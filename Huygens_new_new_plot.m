clc
clear all


N = 3;                %number of Fourier orders
L = 2;                 %number of layers
periodx = 360;  %period of periodic layer
periody = 360;  %period of periodic layer
r = 130;        %disc radius
a = periodx;  
h = zeros(L,1);
h(2) = 250;
h(1) = 360;       %thickness of periodic layer

M = 501;               %number of modes for Fourier transform of epsilon
Mr = (r/a)*M;

i0 = 1+floor(M/2);
j0 = 1+floor(M/2);
%{
a = 360
H = 360
R = 130
n1 = n2 = 1.46

???????? ?? 1300 ?? 1700, ???? ?? 35 ?? 80 ????????.

????? 4 ???????? ? ??????? eps: ????????? ? ???? ??? ?? ? ?? ???????????.

%}
lmin = 1300;
lmax = 1700;
Nl=81;
lambda = linspace(lmin,lmax,Nl);

n_media = 1.46;
eps_media = n_media^2;
n_prism = 2.3;
eps_prism = n_prism^2;
n_substrate = 1.46;

Si_dispersion = xlsread('silicon_cryst_500-1500nm.xlsx');
Si_lambda = Si_dispersion(:,1)*1000;

n_Si = zeros(Nl,1);
eps_Si = zeros(Nl,1);


for i=1:Nl
    [ll,num] = min( abs (lambda(i)-Si_lambda(:) ) );
    %llambda = lambda(i)
    %si_llambda = Si_lambda(num)
    n_Si(i) = Si_dispersion(num,2)+ 1j*Si_dispersion(num,3);
    eps_Si(i) = Si_dispersion(num,5) + 1j*Si_dispersion(num,6);
end




thetamin = 35*pi/180;
thetamax = 80*pi/180;
Nt=91;
theta = linspace(thetamin,thetamax,Nt);
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


[num1,num2] = size(Si_dispersion);
eps_Si_final = Si_dispersion(num1,5) + 1j*Si_dispersion(num1,6);

epsilon = zeros(M, M, L);
refIndices = [n_prism n_substrate];
epsilon(:,:,2) = eps_media*ones(M,M);
for nlayer=2:2
    [eps11(:,:,nlayer), eps22(:,:,nlayer), eps33(:,:,nlayer)] =...
        FMM_eps123_new(epsilon(:,:,nlayer),N,M);
end
%{
for ii=1:M
    for jj=1:M
        if ( ((ii-i0)^2+(jj-j0)^2) <= Mr^2)
            epsilon(jj,ii,1) = 12.2+1j*0.5;
        else
            epsilon(jj,ii,1) = eps_media;
            
        end
    end
end
for nlayer=1:L
    [eps11(:,:,nlayer), eps22(:,:,nlayer), eps33(:,:,nlayer)] =...
        FMM_eps123_new(epsilon(:,:,nlayer),N,M);
end
%}
%{
for ii=1:M
    for jj=1:M
        if ( ((ii-i0)^2+(jj-j0)^2) <= Mr^2)
            epsilon(jj,ii,1) = eps_Si_final; %eps_Si(i);
        else
            epsilon(jj,ii,1) = eps_media;
            
        end
    end
end
for nlayer=1:L
    [eps11(:,:,nlayer), eps22(:,:,nlayer), eps33(:,:,nlayer)] =...
        FMM_eps123_new(epsilon(:,:,nlayer),N,M);
end
%}

phase_R = zeros(Nl,Nt);
Rsum_full = zeros(Nl,Nt);
polarization = 'TM';
load('TE_a_360_H_360_R_130_n1_1_46_n2_1_46.mat', 'Rsum', 'phase_R', 'lambda', 'theta')

%{
a = 360
H = 360
R = 130
n1 = n2 = 1.46

???????? ?? 1300 ?? 1700, ???? ?? 35 ?? 80 ????????.
%}
%save('MyMatrix.txt', 'A', '-ascii', '-double', '-tabs')
%llambda = transpose(lambda);
%ttheta = transpose(theta);
%lt = cat(2,llambda,ttheta);
%{
save('TM_R_a_360_H_360_R_130_n1_1_46_n2_1_46.txt', 'Rsum', '-ascii', '-double', '-tabs')
save('TM_phase_a_360_H_360_R_130_n1_1_46_n2_1_46.txt', 'phase_R','-ascii',  '-double', '-tabs')
save('TM_lambda_theta_a_360_H_360_R_130_n1_1_46_n2_1_46.txt', 'lambda', 'theta','-ascii', '-double', '-tabs')
save('TM_a_360_H_360_R_130_n1_1_46_n2_1_46.mat', 'Rsum', 'phase_R', 'lambda', 'theta')
%}
%{
figure(1)
plot(lambda, Rsum, 'g', lambda, Tsum, 'r', 'Linewidth', 2);
%plot(lambda, transpose(Rsum), 'r', 'Linewidth', 2);
hold off
%}
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



c = physconst('LightSpeed');
h = 4.135666 * 10^(-15);
kx = (2*a*n_prism./lambda)'*sin(theta);
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

plot(frequency*n_prism*a*2*10^3/c,frequency,'b',...
    frequency*n_substrate*a*2*10^3/c,frequency,'k','Linewidth',4)
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
caxis([-0.1 0.1])
hcb=colorbar
title(hcb,'phase, deg')
hold on

plot(frequency*n_prism*a*2*10^3/c,frequency,'b',...
    frequency*n_substrate*a*2*10^3/c,frequency,'k','Linewidth',4)
hold off
