clc
clear all


N = 3;                %number of Fourier orders
L = 2;                 %number of layers
periodx = 300;  %period of periodic layer
periody = 300;  %period of periodic layer
r = 100;        %disc radius
a = periodx;  
h = zeros(L,1);
h(2) = 300;
h(1) = 100;       %thickness of periodic layer

M = 701;               %number of modes for Fourier transform of epsilon
Mr = (r/a)*M;


i0 = 1+floor(M/2);
j0 = 1+floor(M/2);

lmin = 1000;
lmax = 1800;
Nl=40;
lambda = linspace(lmin,lmax,Nl);

n_media = 1.0;%1.66;
eps_media = n_media^2;
n_substrate = 2.13;
n_prism = 2.4;
%eps_prism = n_prism^2;

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

epsilon = zeros(M, M, L, Nl);
for i_lambda = 1:Nl
    epsilon(:,:,2,i_lambda) = eps_media*ones(M,M);
end
for i=1:M
    for j=1:M
        if ( ((i-i0)^2+(j-j0)^2) <= Mr^2)
            epsilon(j,i,1,:) = eps_Si(:);
        else
            epsilon(j,i,1,:) = eps_media;
            
        end
    end
end

refIndices = [n_prism n_substrate];


thetamin = 20*pi/180;
thetamax = 60*pi/180;
Nt=40;
theta = linspace(thetamin,thetamax,Nt);
phi = 0*pi/180;
Np=1;
Rsum=zeros(Nl,Nt);
Tsum=zeros(Nl,Nt);


P = 2*N+1;
Q = 2*N+1;
R = 1;

eps11=zeros(P*Q,P*Q,L,Nl);
eps22=zeros(P*Q,P*Q,L,Nl);
eps33=zeros(P*Q,P*Q,L,Nl);
for i=1:L
    for k=1:Nl
    [eps11(:,:,i,k), eps22(:,:,i,k), eps33(:,:,i,k)] = FMM_eps123_new(epsilon(:,:,i,k),N,M);
    end
end

polarization = 'TE';
for i=1:Nl
    for j=1:Nt
    for k=1:Np 
        eps11_t(:,:,:) = eps11(:,:,:,i);
        eps22_t(:,:,:) = eps22(:,:,:,i);
        eps33_t(:,:,:) = eps33(:,:,:,i);
    [eta_R, eta_T] = FMM_1D_TE_RT_multi(eps11_t,eps22_t,eps33_t,...
        periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L, polarization);
    lambda(i)
    Rsum(i,j) = sum(eta_R);
    Tsum(i,j) = sum(eta_T);
    end
    end
end

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
xlabel('kx, \pi/a');
ylabel('frequency, THz');
colormap('jet');
colorbar;
set(gca,'fontsize', 16)
shading flat
caxis([0 1])
colorbar
hold on


plot(frequency*n_prism*a*2*10^3/c,frequency,'b',frequency*n_media*a*2*10^3/c,frequency,'k','Linewidth',4)
hold off
