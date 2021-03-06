clc
clear all

N = 3;                %number of Fourier orders
L_real = 5;
L = L_real + (L_real-1); %number of layers considering intersections



periodx = 1300*10^(-9);
periody = 1000*10^(-9);

w_wide = 1050*10^(-9);
w_narrow = 350*10^(-9);

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


M = 2600;

n_SiO2 = 1.45;
epsSiO2 = n_SiO2^2;

n_resist = 1.53;
eps_resist = n_resist^2;

n_air = 1.0;
eps_air = 1.0;

layers(L,:) = 'waveguide';
layers(L-1,:) = 'top_cross';
for i = 1:(L-2)
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
epsilon_top_cross = eps_air*ones(M,M);
refIndices = [n_air n_SiO2];

%in epsilon row is y, column is x
%lowest layer has nlayer=1, top layer has nlayer=L
x = (1:1:M)*periodx/M;
y = (1:1:M)*periody/M;

for i=1:M
    for j=1:M
        %waveguide
        if (x(i)<=w_wide)
            epsilon_waveguide(j,i) = eps_resist;
        end
        %top_cross
        if (x(i)<=w_wide)
            epsilon_top_cross(j,i) = eps_resist;
        end
        if (y(j)<=w_narrow)
            epsilon_top_cross(j,i) = eps_resist;
        end
        %yperiodic
        if (y(j)<=w_narrow)
            epsilon_yperiodic(j,i) = eps_resist;
        end
        %xperiodic
        if (x(i)>(w_wide/2-w_narrow/2) ) && (x(i)<=(w_wide/2+w_narrow/2) )
            epsilon_xperiodic(j,i) = eps_resist;
        end
        %xycrossed
        if (x(i)>(w_wide/2-w_narrow/2) ) && (x(i)<=(w_wide/2+w_narrow/2) )
            epsilon_xycrossed(j,i) = eps_resist;
        end
        if (y(j)<=w_narrow)
            epsilon_xycrossed(j,i) = eps_resist;
        end  
    end
end
%{
lmin = 1333*10^(-9);
lmax = 1336*10^(-9);
thetamin = 0.25*pi/180;
thetamax = 1.75*pi/180;
lambda = linspace(lmin, lmax, 200);
[Nll,Nl] = size(lambda);
theta = linspace(thetamin,thetamax,7);
[Ntt,Nt]=size(theta);
%}
lmin = 1346*10^(-9);
lmax = 1346.6*10^(-9);
thetamin = 0.5*pi/180;
thetamax = 1.5*pi/180;
lambda = linspace(lmin, lmax, 100);
[Nll,Nl] = size(lambda);
theta = linspace(thetamin,thetamax,15);
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
[eps11_top_cross, eps22_top_cross, eps33_top_cross] = FMM_eps123_new(epsilon_top_cross,N,M);
%%%%%%%%
optional_top_layer = 'yes';
epsilon_top = eps_resist*ones(M,M);
[eps11_top,eps22_top,eps33_top] = FMM_eps123_new(epsilon_top,N,M);
%%%%%%%%
eps11(:,:,L) = eps11_waveguide;
eps22(:,:,L) = eps22_waveguide;
eps33(:,:,L) = eps33_waveguide;

eps11(:,:,L-1) = eps11_top_cross;
eps22(:,:,L-1) = eps22_top_cross;
eps33(:,:,L-1) = eps33_top_cross;
if L>1
for i=1:L-2
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
end

Rsum_another_top = zeros(Nl,Nt);
for i=1:Nl
    for j=1:Nt
        for k=1:Np
            [eta_R, eta_R_another_top] = FMM_2D_RT_multi(eps11,eps22,eps33,...
                periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L,polarization,...
                layers, optional_top_layer, eps11_top, eps22_top, eps33_top);
            %Rsum(j,k) = sum(eta_R1);
            %Tsum(j,k) = sum(eta_T1);
            Rsum(i,j) = sum(eta_R);
            Rsum_another_top(i,j) = sum(eta_R_another_top);
            
        end
    end
    llambda = lambda(i)
end


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

%{
for i=1:Nt
    f=figure;
    %plot(lambda, Rsum(:,i), 'b', lambda, Rsum_another_top(:,i)+0.5, 'm', 'LineWidth', 2)
    plot(lambda, Rsum(:,i), 'b', 'LineWidth', 2)
    axis tight
    ax = gca;
    ax.XAxis.MinorTick = 'on';
    xlabel('lambda, mkm')
    ylabel('R')
    set(gca,'fontsize', 16)
end
%}
llambda=transpose(lambda);
data = cat(2,llambda,Rsum(:,1));

for i=2:Nt
    data = cat(2,data,Rsum(:,i));
end
%save('5_layers_05_to_15.mat','data');




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

tl = linspace(lmin,lmax,400);
tt = linspace(thetamin, thetamax,30)*180/pi;
[XI,YI] = meshgrid(tl,tt);
ZI = griddata(theta*180/pi,lambda,Rsum_another_top,YI,XI);
figure;
pcolor(XI,YI,ZI)
shading flat
colormap('jet');
colorbar;
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
%pcolor(theta*180/pi,lambda,Rsum);
%imagesc(theta*180/pi,lambda,Rsum);
%set(gca,'Yscale','linear','Ydir','normal');

