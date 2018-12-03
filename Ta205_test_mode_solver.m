clc
clear all

N = 5;                %number of Fourier orders
periodx = 490*10^(-9);
dx = 145*10^(-9);
periody = 700*10^(-9);
L = 2;
h = zeros(L,1);
h(2) = 90*10^(-9);
h(1) = 1*10^(-6);
M = 5001;
x = (1:1:M)*periodx/M;
epsilon = zeros(M, M, L);


nTa2O5 = 2.1;
epsTa2O5 = nTa2O5^2;
nSiO2 = 1.48;
epsSiO2 = nSiO2^2;
nSi = 3.71;
epsilon(:,:,2)=epsTa2O5*ones(M,M,1);
epsilon(:,:,1)=epsSiO2*ones(M,M,1);
refIndices = [1.0 nSi];  
for i=1:M
    for j=1:M
    if x(i)<dx
        epsilon(j,i,2) = 1.0;
    end
    end
end
%{
lmin = 735*10^(-9);
lmax = 765*10^(-9);
lambda = linspace(lmin, lmax, 80);
[Nll,Nl] = size(lambda);
thetamin = 0*pi/180;
thetamax = 1*pi/180;
theta = linspace(thetamin,thetamax,20);
[Ntt,Nt] = size(theta);
%}
l1 = 760*10^(-9);
dl = 0.05*10^(-9);
l2 = l1+dl;
lambda = [l1 l2];
c = 3*10^8;
w1 = c/l1;
w2 = c/l2;
dw = w2-w1;
Nl = 2;

theta = 4*pi/180;
Nt = 1;

phi = 0*pi/180;
Np=1;


P = 2*N+1;
Q = 2*N+1;
R = 1;

eps11=zeros(P*Q,P*Q,L);
eps22=zeros(P*Q,P*Q,L);
eps33=zeros(P*Q,P*Q,L);
for i=1:L
[eps11(:,:,i), eps22(:,:,i), eps33(:,:,i)] = FMM_eps123_new(epsilon(:,:,i),N,M);
end
polarization = 'TE';

N_iterations = 10;

for ii=1:N_iterations
    
    for i=1:Nl
        for j=1:Nt
            for k=1:Np
                [R] = FMM_1D_TE_RT_multi_mode_solver(eps11,eps22,eps33,...
                    periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L, polarization);
            
                if i==1
                    R0 = R;
                else
                    R1 = R;
                end
            end
        end
    end
    dR = (R1-R0)/dw;
    MAXDR(ii) = max(dR(:));
            
    A=-dR\R0;
    [V, D] = eig(A);
            
    w_array = diag(D);
    dw0 = min(w_array);
            
    w1 = w1 + dw0;
    w2 = w1 + dw;
            
    l1 = c/w1;
    l2 = c/w2;
    lambda = [l1 l2];
            
    dw00(ii) = dw0;
end
l1

