clc
clear all

N = 4;                %number of Fourier orders

L = 1;
h = zeros(L,1);
th = 212*10^(-9);
h(1) = th;
periodx = 1.8*th;
dx = periodx*0.45;
periody = 1.8*th;
M = 3001;
x = (1:1:M)*periodx/M;
epsilon = zeros(M, M, L);


nlattice = 2.08;
epslattice = nlattice^2;
nmedia = 1.46;
epsmedia = nmedia^2;
epsilon(:,:,1)=epsmedia*ones(M,M,1);
refIndices = [nmedia nmedia];  
for i=1:M
    for j=1:M
    if x(i)<dx
        epsilon(j,i,1) = epslattice;
    end
    end
end

l1 = 760*10^(-9);
ll1 = l1;
dl = (0.05+0.1j)*10^(-9);
l2 = l1+dl;
lambda = [l1 l2];
c = 3*10^8;
w1 = c/l1;
w2 = c/l2;
dw = w2-w1;
Nl = 2;

kx = 9.06*10^6;
%theta = asin(kx*l1/(2*pi*nmedia));
theta1 = 18*pi/180;
thetamin = theta1;% - 2*pi/180;
thetamax = theta1 + 0.5*pi/180;
%Nt = 1;
theta = linspace(thetamin,thetamax,10);
[Ntt,Nt] = size(theta);


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
N_iterations = 5;

w_eig = zeros(Nt,1);
dw00 = zeros(N_iterations, Nt);
dw0 = w1;

for j=1:Nt
    j
    for k=1:Np
        %ii=0
        for ii=1:N_iterations
        %while ((abs(dw0/w1)>10^(-5)) && ii<5)
        %    ii=ii+1
        ii
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
            
            w_array = sort(diag(D))
            dw0 = min(w_array);
            if ii==1
                dw0 = w_array(31);
            end
            w1 = w1 + dw0;
            w2 = w1 + dw;
            
            l1 = c/w1;
            l2 = c/w2;
            lambda = [l1 l2];
            
            dw00(ii,j) = dw0;
            
        end
        w_eig(j) = w1;
        lambda = [ll1 ll1+dl];
    end
    
end
lambda_eig = c./w_eig;
k0=2*pi./real(lambda_eig);
%kx = k0*periodx*sin(theta')*cos(phi);
%plot(kx,w_eig,'b','Linewidth',2);

theta_new = theta*180/pi;
%plot(theta_new, lambda_eig*10^9)
plot(lambda_eig*10^9,theta_new)
theta_new1 = theta_new(2:10);
lambda_eig1 = lambda_eig(2:10)*10^9;
plot(lambda_eig1*10^9,theta_new1)

p=polyfit(transpose(theta_new1),lambda_eig1,4);
x1 = theta_new1;
y1=polyval(p,x1);
figure(1)
plot(lambda_eig1,theta_new1,'o')
hold on
plot(y1,x1,'r')
hold off
%{
k0=2*pi/lambda;
kx = sin(theta')*cos(phi);
ky = sin(theta')*sin(phi);
figure;
pcolor(kx,ky,Rsum)
%shading interp
colormap gray;
%}
