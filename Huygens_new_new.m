clc
clear all


N = 3;                %number of Fourier orders
L = 2;                 %number of layers
periodx = 300;  %period of periodic layer
periody = 300;  %period of periodic layer
r = 125;        %disc radius
a = periodx;  
h = zeros(L,1);
h(2) = 300;
h(1) = 375;       %thickness of periodic layer

M = 501;               %number of modes for Fourier transform of epsilon
Mr = (r/a)*M;

i0 = 1+floor(M/2);
j0 = 1+floor(M/2);

%%%%%%%%%%%%%%%%%%%%
%{
a = 300 nm
H = 375 nm
R = 125 nm
n1 = n2 = 1.46

Диапазон длин волн от 1100 до 1600 нм (50 шагов), угол от 25 до 60 градусов (35 шагов).
%}
lmin = 600;
lmax = 1100;
Nl=51;
lambda = linspace(lmin,lmax,Nl);

n_media = 1.46;
eps_media = n_media^2;
n_prism = 2.5;
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




thetamin = 25*pi/180;
thetamax = 60*pi/180;
Nt=36;
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
polarization = 'TE';

for i=1:Nl
    
    for ii=1:M
        for jj=1:M
            if ( ((ii-i0)^2+(jj-j0)^2) <= Mr^2)
                epsilon(jj,ii,1) = eps_Si(i);
            else
                epsilon(jj,ii,1) = eps_media;
                
            end
        end
    end
    for nlayer=1:L
        [eps11(:,:,nlayer), eps22(:,:,nlayer), eps33(:,:,nlayer)] =...
            FMM_eps123_new(epsilon(:,:,nlayer),N,M);
    end
    
    for j=1:Nt
        for k=1:Np
            
            [eta_R, eta_T] = FMM_1D_TE_RT_multi(eps11,eps22,eps33,...
                periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L, polarization);
            
            Rsum(i,j) = sum(eta_R);
            Tsum(i,j) = sum(eta_T);
        end
    end
    lambda(i)
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


plot(frequency*n_prism*a*2*10^3/c,frequency,'b',frequency*n_substrate*a*2*10^3/c,frequency,'k','Linewidth',4)
hold off
