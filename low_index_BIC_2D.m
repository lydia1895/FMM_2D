clc
clear all

N = 3;                %number of Fourier orders
L = 10; %number of layers

periodx = 1300*10^(-9);
periody = 1000*10^(-9);

w_wide = 1050*10^(-9);
w_narrow = 350*10^(-9);
w=w_narrow*ones(L,1);
w(L) = w_wide;

h = zeros(L,1);
h_standard = 800*10^(-9);
h = h_standard*ones(L,1);
%h(1) = 2*h_standard;
%h(2) = 2*h_standard;


M = 2601;

epsilon = zeros(M, M, L);


n_SiO2 = 1.45;
epsSiO2 = n_SiO2^2;

n_resist = 1.53;
eps_resist = n_resist^2;

n_air = 1.0;
eps_air = 1.0;

layers(L,:) = 'waveguide';
for i = 1:(L-1)
    if mod(L-1-i,2) == 0
        layers(i,:) = 'yperiodic';
    else
        layers(i,:) = 'xperiodic';
    end
end
        

epsilon = eps_air*ones(M,M,L);
refIndices = [n_air n_SiO2];

%in epsilon row is y, column is x
%lowest layer has nlayer=1, top layer has nlayer=L
x = zeros(M,L);
y = zeros(M,L);
for i=1:L
    x(:,i) = (1:1:M)*periodx/M;
    y(:,i) = (1:1:M)*periody/M;
end

for i=1:M
    for j=1:M
        if (x(i,L)<w(L)) %|| (periodx_top_real<x(i,3))&&(x(i,3)<=periodx_top_real+wx(3)) %||...
                %(2*periodx_top_real<x(i,3))&&(x(i,3)<2*periodx_top_real+wx(3))
            epsilon(j,i,L) = eps_resist;
        end
        if y(j,L-1)<w(L-1)
            epsilon(j,i,L-1) = eps_resist;
        end
        if (x(i,L-2)>(w(L)/2-w(L-2)/2) ) && (x(i,L-2)<(w(L)/2+w(L-2)/2) )
            epsilon(j,i,L-2) = eps_resist;
        end
    end
end
for i = 1:(L-3)
    if mod(abs(i-(L-1)),2) == 0
        epsilon(:,:,i)=epsilon(:,:,L-1);
    else
        epsilon(:,:,i)=epsilon(:,:,L-2);
    end
end

lmin = 1410*10^(-9);
lmax = 1510*10^(-9);
thetamin = 1.5*pi/180;
thetamax = 5*pi/180;
lambda = linspace(lmin, lmax, 300);
[Nll,Nl] = size(lambda);
theta = linspace(thetamin,thetamax,8);
[Ntt,Nt]=size(theta);
%theta = [0.1 0.5 1]*pi/180;
%Nt=3;
%theta = 0.3*pi/180;
%Nt=1;

phi = 0*pi/180;
Np=1;
Rsum=zeros(Nl,Nt);
Tsum=zeros(Nl,Nt);

polarization = 'TE';

P = 2*N+1;
Q = 2*N+1;
R = 1;

eps11=zeros(P*Q,P*Q,L);
eps22=zeros(P*Q,P*Q,L);
eps33=zeros(P*Q,P*Q,L);
for i=1:L
    [eps11(:,:,i), eps22(:,:,i), eps33(:,:,i)] = FMM_eps123_new(epsilon(:,:,i),N,M);
end



for i=1:Nl
    for j=1:Nt
        for k=1:Np
            [eta_R1, eta_T1] = FMM_2D_RT_multi(eps11,eps22,eps33,...
                periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L,polarization, layers);
            %Rsum(j,k) = sum(eta_R1);
            %Tsum(j,k) = sum(eta_T1);
            Rsum(i,j) = sum(eta_R1);
            Tsum(i,j) = sum(eta_T1);
            
        end
    end
    llambda = lambda(i)
end
%{
epsilon(:,:,2)=epsTa2O5*ones(M,M,1);

Rsum_non=zeros(Nl,Nt);
Tsum_non=zeros(Nl,Nt);

eps11=zeros(P*Q,P*Q,L);
eps22=zeros(P*Q,P*Q,L);
eps33=zeros(P*Q,P*Q,L);
for i=1:L
[eps11(:,:,i), eps22(:,:,i), eps33(:,:,i)] = FMM_eps123_new(epsilon(:,:,i),N,M);
end



for i=1:Nl
    for j=1:Nt
    for k=1:Np
    [eta_R1, eta_T1] = FMM_2D_RT_multi(eps11,eps22,eps33,...
        periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L,polarization);
    %Rsum(j,k) = sum(eta_R1);
    %Tsum(j,k) = sum(eta_T1);
    Rsum_non(i,j) = sum(eta_R1);
    Tsum_non(i,j) = sum(eta_T1);
    end
    end
end
%}

lambda = lambda*10^6;
lmin = lmin*10^6;
lmax = lmax*10^6;



figure(1)
%{
plot(lambda, Rsum(:,1), 'b', lambda, Rsum(:,2)+0.5, 'g',...
    lambda, Rsum(:,3)+1.0,'r','LineWidth', 2)
    %lambda, Rsum(:,4)+1.5, 'm', 'LineWidth', 2)
h5 = legend('theta=0.1 deg','theta=0.5 deg','theta=1 deg',3);
%}
for i=1:Nt
    f=figure;
    plot(lambda, Rsum(:,i), 'b', 'LineWidth', 2)
    %set(h5,'Interpreter','none')
    axis tight
    ax = gca;
    ax.XAxis.MinorTick = 'on';
    xlabel('lambda, mkm')
    ylabel('R')
    set(gca,'fontsize', 16)
end

llambda=transpose(lambda);
data = cat(2,llambda,Rsum(:,1));

for i=2:Nt
    data = cat(2,data,Rsum(:,i));
end
save('10_layers_1.5_to_3_deg.mat','data');


%{
figure(2)
tl = linspace(lmin,lmax,400);
tt = linspace(thetamin, thetamax,40)*180/pi;
[XI,YI] = meshgrid(tl,tt);
ZI = griddata(theta*180/pi,lambda,Rsum,YI,XI);
figure;
%pcolor(kx,ky,Rsum)
pcolor(XI,YI,ZI)
shading flat
caxis([0 1])
colormap('jet');
colorbar;
hold off
%}


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

