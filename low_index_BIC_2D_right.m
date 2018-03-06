clc
clear all

N = 3;                %number of Fourier orders
L_real = 11;
L = L_real + (L_real-1); %number of layers considering intersections

periodx = 1300*10^(-9);
periody = 1000*10^(-9);

w_wide = 1050*10^(-9);
w_narrow = 350*10^(-9);
w=w_narrow*ones(L,1);
w(L) = w_wide;


h_layer = 600*10^(-9);
h_intersection = 100*10^(-9);
h = h_layer*ones(L,1);
h(1) = 700*10^(-9);
h(L) = 700*10^(-9);
for i=2:L-1
    if mod(i,2)==0
        h(i) = h_intersection;
    end
end


M = 5200;

n_SiO2 = 1.45;
epsSiO2 = n_SiO2^2;

n_resist = 1.53;
eps_resist = n_resist^2;

n_air = 1.0;
eps_air = 1.0;

layers(L,:) = 'waveguide';
for i = 1:(L-1)
    if mod(i,2) == 1
        i_real = (i+1)/2;
        if mod(L_real-1-i_real,2)==0
            layers(i,:) = 'yperiodic';
        else
            layers(i,:) = 'xperiodic';
        end
    else
        layers(i,:) = 'xycrossed';
    end
end
        

epsilon_waveguide = eps_air*ones(M,M);
epsilon_yperiodic = eps_air*ones(M,M);
epsilon_xperiodic = eps_air*ones(M,M);
epsilon_xycrossed = eps_air*ones(M,M);
refIndices = [n_air n_SiO2];

%in epsilon row is y, column is x
%lowest layer has nlayer=1, top layer has nlayer=L
x = (1:1:M)*periodx/M;
y = (1:1:M)*periody/M;

for i=1:M
    for j=1:M
        %waveguide
        if (x(i)<=w(L))
            epsilon_waveguide(j,i) = eps_resist;
        end
        %yperiodic
        if (y(j)<=w(L-1))
            epsilon_yperiodic(j,i) = eps_resist;
        end 
        %xperiodic
        if (x(i)>(w(L)/2-w(L-2)/2) ) && (x(i)<=(w(L)/2+w(L-2)/2) )
            epsilon_xperiodic(j,i) = eps_resist;
        end
        %xycrossed    
        if (x(i)>(w(L)/2-w(L-2)/2) ) && (x(i)<=(w(L)/2+w(L-2)/2) )
            epsilon_xycrossed(j,i) = eps_resist;
        end
        if (y(j)<=w(L-1))
            epsilon_xycrossed(j,i) = eps_resist;
        end        
    end
end

lmin = 1465*10^(-9);
lmax = 1480*10^(-9);
thetamin = 0.05*pi/180;
thetamax = 0.50*pi/180;
lambda = linspace(lmin, lmax, 300);
[Nll,Nl] = size(lambda);
theta = linspace(thetamin,thetamax,10);
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
[eps11_waveguide, eps22_waveguide, eps33_waveguide] = FMM_eps123_new(epsilon_waveguide,N,M);
[eps11_yperiodic, eps22_yperiodic, eps33_yperiodic] = FMM_eps123_new(epsilon_yperiodic,N,M);
[eps11_xperiodic, eps22_xperiodic, eps33_xperiodic] = FMM_eps123_new(epsilon_xperiodic,N,M);
[eps11_xycrossed, eps22_xycrossed, eps33_xycrossed] = FMM_eps123_new(epsilon_xycrossed,N,M);

eps11(:,:,L) = eps11_waveguide;
eps22(:,:,L) = eps22_waveguide;
eps33(:,:,L) = eps33_waveguide;
for i=1:L-1
    if strcmp(layers(i,:),'yperiodic')==1
        eps11(:,:,i) = eps11_yperiodic;
        eps22(:,:,i) = eps22_yperiodic;
        eps33(:,:,i) = eps33_yperiodic;
    end
    if strcmp(layers(i,:),'xperiodic')==1
        eps11(:,:,i) = eps11_xperiodic;
        eps22(:,:,i) = eps22_xperiodic;
        eps33(:,:,i) = eps33_xperiodic;
    end
    if strcmp(layers(i,:),'xycrossed')==1
        eps11(:,:,i) = eps11_xycrossed;
        eps22(:,:,i) = eps22_xycrossed;
        eps33(:,:,i) = eps33_xycrossed;
    end  
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
%{
figure(2)
hold on
plot(lambda, Rsum(:,1), 'b', lambda, Rsum(:,2)+0.5, 'g',...
    lambda, Rsum(:,3)+1, 'r', 'LineWidth', 2)
h5 = legend('theta=0.1','theta=0.5','theta=1',3);
set(h5,'Interpreter','none')
%title('W=',w*10^9,' nm, D=', period*10^9,' nm')
axis tight
axis([lmin lmax 0 2.5])
set(gca,'fontsize', 18)
hold off
%}
%{
plot(lambda, Rsum(:,1), 'b', lambda, Rsum(:,2)+0.5, 'g',...
    lambda, Rsum(:,3)+1.0,'r','LineWidth', 2)
    %lambda, Rsum(:,4)+1.5, 'm', 'LineWidth', 2)
h5 = legend('theta=0.1 deg','theta=0.5 deg','theta=1 deg',3);
%}


for i=1:Nt
    f=figure;
    plot(lambda, Rsum(:,i), 'b', 'LineWidth', 2)
    set(h5,'Interpreter','none')
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
save('11_layers_right.mat','data');



%{
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

