clc
clear all

N=3;
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


lmin = 330*10^(-9);
lmax = 370*10^(-9);
lambda = linspace(lmin, lmax, 70);
[Nll,Nl] = size(lambda);


kx = 9.06*10^6;
l1 = 489*10^(-9);
theta1 = asin(kx*l1/(2*pi))
thetamin = 0;
thetamax = 5*pi/180;
theta = linspace(thetamin,thetamax,20);
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
plot(lambda, Rsum(:,7), 'b', 'LineWidth', 2)
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

