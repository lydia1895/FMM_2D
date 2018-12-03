clc
clear all

N = 3;                %number of Fourier orders
periodx = 490*10^(-9);
dx = 145*10^(-9);
periody = 650*10^(-9);
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

lmin = 690*10^(-9);
lmax = 760*10^(-9);
lambda = linspace(lmin, lmax, 30);
[Nll,Nl] = size(lambda);

thetamin = 0*pi/180;
thetamax = 5*pi/180;
theta = linspace(thetamin,thetamax,15);
[Ntt,Nt] = size(theta);

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
for i=1:L
[eps11(:,:,i), eps22(:,:,i), eps33(:,:,i)] = FMM_eps123_new(epsilon(:,:,i),N,M);
end
polarization = 'TE';
for i=1:Nl
    for j=1:Nt
    for k=1:Np 
    [eta_R1, eta_T1] = FMM_1D_TE_RT_multi(eps11,eps22,eps33,...
        periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L, polarization);
    %Rsum(j,k) = sum(eta_R1);
    %Tsum(j,k) = sum(eta_T1);
    Rsum(i,j) = sum(eta_R1);
    Tsum(i,j) = sum(eta_T1);
    end
    end
end


lambda = lambda*10^6;
lmin = lmin*10^6;
lmax = lmax*10^6;

figure(5)
hold on
plot(lambda, Rsum(:,1), 'b', 'LineWidth', 2)
%plot(lambda, Rsum(:,1), 'b',  lambda, Rsum(:,2)+1, 'r', 'LineWidth', 2)

%h5 = legend('theta=0','theta=1','theta=2','theta=3',4);

%title('W=',w*10^9,' nm, D=', period*10^9,' nm')
axis tight
%axis([lmin lmax 0 2.5])
set(gca,'fontsize', 18)
xlabel('lambda, mkm')
ylabel('R')
%h5 = legend('theta=0','theta=0.2','theta=0.4','theta=0.6',4);
%set(h5,'Interpreter','none')
hold off
%plot(lambda, Rsum)
%{
Collist = 'ybmgr'
col = randi([1 5])
figure(2)
hold on
plot(lambda, Rsum+0.5, Collist(col))
hold off
%}
%plot(lambda, Rsum)
%plot(lambda, Rsum(:,1,1), 'b', lambda, Rsum(:,2,1), 'g', lambda, Rsum(:,3,1), 'r');

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
tl = linspace(lmin,lmax,400);
tt = linspace(thetamin, thetamax,30)*180/pi;
[XI,YI] = meshgrid(tl,tt);
ZI = griddata(theta*180/pi,lambda,Rsum,YI,XI);
figure;
pcolor(XI,YI,ZI)
shading flat
colormap('jet');
colorbar;
hold off

%plot([0:1:TH-1],Rsum,'b',[0:1:TH-1],Tsum,'r')

%num = zeros(2*N+1,1);
%for q=1:(2*N+1)
%   num(q) = q-N-1;
%end
%bar(num,eta,'stack')
%bar(num, eta_T1)