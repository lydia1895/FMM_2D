clc
clear all

N = 6;                %number of Fourier orders

L = 1;
h = zeros(L,1);
th = 212*10^(-9);
h(1) = th;
periodx = th;
dx = periodx*0.45;
periody = th;
M = 3001;
x = (1:1:M)*periodx/M;
epsilon = zeros(M, M, L);


nlattice = 2.125;
epslattice = nlattice^2;
nmedia = 1.46;
epsmedia = nmedia^2;
epsilon(:,:,1)=epsmedia*ones(M,M,1);
refIndices = [nmedia nmedia];  
for i=1:M
    for j=1:M
    if x(i)<dx
        epsilon(j,i,2) = epslattice;
    end
    end
end

l1 = 489*10^(-9);
dl = 0.05*10^(-9);
l2 = l1+dl;
lambda = [l1 l2];
c = 3*10^8;
w1 = c/l1;
w2 = c/l2;
dw = w2-w1;
Nl = 2;

kx = 9.06*10^6;
theta = asin(kx*l1/(2*pi));
%thetamax = thetamin + 1*pi/180;
Nt = 1;
%theta = linspace(thetamax,thetamin,6);
%[Ntt,Nt] = size(theta);


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
polarization = 'TM';
N_iterations = 12;

w_eig = zeros(Nt,1);
dw00 = zeros(N_iterations, Nt);

for j=1:Nt
    for k=1:Np
        for ii=1:N_iterations
    
            for i=1:Nl
        
                [R] = FMM_1D_TE_RT_multi_mode_solver(eps11,eps22,eps33,...
                    periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L, polarization);
            
                if i==1
                    R0 = R;
                else
                    R1 = R;
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
            
            dw00(ii,j) = dw0;
        end
        w_eig(j) = w1;
    end
    
end
lambda_eig = c./w_eig;
k0=2*pi./real(lambda_eig);
kx = k0*periodx*sin(theta')*cos(phi);
plot(kx,w_eig,'b','Linewidth',2);

theta_new = theta*180/pi;
%plot(theta_new, lambda_eig*10^9)
plot(lambda_eig*10^9,theta_new)
%{
k0=2*pi/lambda;
kx = sin(theta')*cos(phi);
ky = sin(theta')*sin(phi);
figure;
pcolor(kx,ky,Rsum)
%shading interp
colormap gray;
%}
