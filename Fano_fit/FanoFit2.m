%function [ x ] = FanoFit2( data, width,q,H,A,B,C)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

<<<<<<< HEAD
% Find all windows of type figure, which have an empty FileName attribute.
allPlots = findall(0, 'Type', 'figure', 'FileName', []);
% Close.
delete(allPlots);

load('11_layers_0.55_to_1.0_new.mat','data')
N=length(data(1,:))
xdata=data(:,1);%xdata=data(:,1);
Q=zeros(N-1,1);
theta = linspace(0.05,0.50,N-1);

qf = linspace(3,4,N-1);
qf = [1100 1000 8.7 8.0 5 5 4 3.6 3.5 3.4];

Hf = linspace(0.3,0.1,N-1);
widthf = linspace(0.001,0.002,N-1);

for i=1:N-1
%i=2;

ydata=data(:,i+1)-min(data(:,i+1));
=======
load('12_layers_BIC.mat','data')
N=length(data(1,:))


thetamin = 0.5;
thetamax = 5;
theta_full = linspace(thetamin,thetamax,10);

xxdata=data(:,1);
left = 1.332;
right = 1.3358;
[elem,num_elem]=min(abs(xxdata-left));
[elem2,num_elem_2]=min(abs(xxdata-right));

ddata = data(num_elem:num_elem_2,:);

i=10;

xdata=ddata(:,1);
ydata=ddata(:,i+1)-min(ddata(:,i+1));

>>>>>>> d5e91902445f43678413472a44e23d3000a77a65
deriv = diff(ydata)./diff(xdata);
    [~,ind]= min(deriv);
    center = xdata(ind);
    
<<<<<<< HEAD
    width = widthf(i);
    q = qf(i);
    H = Hf(i);
=======
    width = 0.0001;
    q = -90;
    H = 0.6;
>>>>>>> d5e91902445f43678413472a44e23d3000a77a65
    A = 0;
    B = 0;
    C = 0;
%       w0    w    asymm
x0 = [center,width,q,H,A,B,C];
%for r/L=0.7 x0 = [center,width,q,0,0.3,2.25];
%Fano = a+bx+c*(q+2(x-x0)/w)^2/(1+4(x-x0)^2/w^2)
Fano = @(x,xdata)(x(7)+x(6)*xdata+x(5)*xdata.^2+...
    x(4)*(x(3)+2*(xdata-x(1))/x(2)).^2./(1+4*(xdata-x(1)).^2/(x(2)^2)));

x=lsqcurvefit(Fano,x0,xdata,ydata);
  f=figure;
  plot(xdata,ydata,'-k',xdata,Fano(x,xdata),'-r');
<<<<<<< HEAD
  Q(i) = abs(x(1)/x(2));
end

theta_2=theta;
Q_2=Q;
save('11_layers_0.55_to_1.0_new.mat','data','theta_2','Q_2')
g=figure;
plot(theta,Q,'-sg','LineWidth',2);
=======
  hold off
  Q(i) = abs(x(1)/x(2));


%g=figure
plot(theta_full,Q,'-sg','Linewidth',2)
theta(i) = theta_full(i);
save('12_layers_BIC_results_1.mat','theta_full','Q','i')
>>>>>>> d5e91902445f43678413472a44e23d3000a77a65
%end

