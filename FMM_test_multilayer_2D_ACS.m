clc
clear all

N = 3;                %number of Fourier orders
periodx = 1000*(10^(-9));  %period of periodic layer
dx = 150*(10^(-9));     %ridge width
periody = 1000*(10^(-9));  %period of periodic layer
dy = 1000*(10^(-9));     
h = zeros(3,1);
h(1) = 3*10^(-6);
h(2) = 350*(10^(-9));       %thickness of periodic layer
h(3) = 1*10^(-6);
M = 1001;              %number of modes for Fourier transform of epsilon
L = 3;                 %number of layers
x = (1:1:M)*periodx/M;
epsilon = zeros(M, M, L);
for i=1:M
    for j=1:M
    epsilon(i,j,1) = 1.44^2;
    if x(i)<dx
        epsilon(i,j,2) = 3.48^2;
    else
        epsilon(i,j,2) = 1.44^2;
    end
    epsilon(i,j,3) = 1.44^2;
    end
end
refIndices = [1.0 3.48];     


%{
M = 1001;
N = 3;                %number of Fourier orders
periodx = 1000*(10^(-9));  %period of periodic layer
dx = 150*(10^(-9));     %ridge width
periody = 1000*(10^(-9));  %period of periodic layer
dy = 1000*(10^(-9));     
epsilon = zeros(M, 1);
L = 1;
h(1) = 0.5*10^(-6);
for i=1:M
    for j = 1:M
        if ((i*periodx)/M<dx) && ((j*periody)/M<dy)
            epsilon(i,j,1) = 9.0;
        else
            epsilon(i,j,1) = 9.0;
        end
    %epsilon(i,2) = 9.0;
    end
end
%}


l1 = 1400*10^(-9);
l2 = 1500*10^(-9);
nl = 50;
lambda = linspace(l1,l2,nl);
%theta = zeros(100,1);
%theta=zeros(90,1);
theta(1) = 2*pi/180;
phi(1) = 0;

%TH = 10;
Rsum=zeros(nl,1);
Tsum=zeros(nl,1);
%theta=zeros(TH,1);
%phi = 0;

for i=1:1:nl
    for j=1:1:1%180
    %theta(i) = (i-1)*pi/180;
    %phi(j) = (j-1)*pi/180;
    [eta_R1, eta_T1] = FMM_1D_TE_RT_multi(epsilon, periodx, periody, h, lambda(i), theta, phi(j), refIndices, N, M, L);
    Rsum(i) = sum(eta_R1);
    Tsum(i) = sum(eta_T1);
    end
end

plot(lambda, Rsum)
%plot([0:1:TH-1],Rsum,'b',[0:1:TH-1],Tsum,'r')
