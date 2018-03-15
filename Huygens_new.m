clc
clear all


N = 3;                %number of Fourier orders
L = 1;                 %number of layers
periodx = 666;  %period of periodic layer
periody = 666;  %period of periodic layer
r = 242;        %disc radius
a = periodx;  
h = zeros(L,2);
h(2) = 
h(1) = 220;       %thickness of periodic layer
%h(1) = 2*10^(-6);
%h(1) = 10^(-3);
M = 1001;               %number of modes for Fourier transform of epsilon
Mr = (r/a)*M;


i0 = 1+floor(M/2);
j0 = 1+floor(M/2);

lmin = 1200;
lmax = 1500;
Nl=30;
lambda = linspace(lmin,lmax,Nl);

n_media=1.66;

Si_dispersion = xlsread('silicon_cryst_500-1500nm.xlsx');
Si_lambda = Si_dispersion(:,1)*1000;
n_Si = zeros(Nl,1);
eps_Si = zeros(Nl,1);


for i=1:Nl
    [ll,num] = min( abs (lambda(i)-Si_lambda(:) ) );
    %llambda = lambda(i)
    %si_llambda = Si_lambda(num)
    n_Si(i) = Si_dispersion(num,2); %+ 1j*Si_dispersion(num,3);
    eps_Si(i) = Si_dispersion(num,5);% + 1j*Si_dispersion(num,6);
end

epsilon = zeros(M, M, Nl);

for i=1:M
    for j=1:M
        for i_lambda = 1:Nl
            if ( ((i-i0)^2+(j-j0)^2) <= Mr^2)
                epsilon(j,i,i_lambda) = eps_Si(i_lambda);
            else
                epsilon(j,i,i_lambda) = n_media^2;
                
            end
        end
    end
end
refIndices = [n_media n_media];


theta = 0*pi/180;
Nt=1;
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
    [eps11(:,:,i,k), eps22(:,:,i,k), eps33(:,:,i,k)] = FMM_eps123_new(epsilon(:,:,k),N,M);
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
    
    Rsum(i,j) = sum(eta_R);
    Tsum(i,j) = sum(eta_T);
    end
    end
end


figure(1)
plot(lambda, Rsum, 'g', lambda, Tsum, 'r', 'Linewidth', 2);
%plot(lambda, transpose(Rsum), 'r', 'Linewidth', 2);
hold off

