%function [ x ] = FanoFit2( data, width,q,H,A,B,C)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Find all windows of type figure, which have an empty FileName attribute.
allPlots = findall(0, 'Type', 'figure', 'FileName', []);
% Close.
delete(allPlots);
clc
clear all

load('11_layers_1.05_to_3.0_new.mat','data')
N=length(data(1,:))
xdata=data(:,1);%xdata=data(:,1);
Q=zeros(N-1,1);
theta = linspace(1.05,3.0,N-1);

qf = [ -5 -5 -5 -5 -5 -5 -6 -10 -90 -100];
qf = linspace(0.9,0,N-1);
Hf = linspace(0.6,0.1,N-1);
Hf(1) = 0.1;
Hf(2) = 0.1;
Hf(3) = 0.2;
qf(1) = 3.4;
widthf = linspace(0.002,0.002,N-1);
  
NN=length(data(:,1));
ddata = data(1:floor(NN/2),:);
xdata=ddata(:,1)    
    
for i=13:N-1
%i=2;

ydata=ddata(:,i+1)-min(ddata(:,i+1));
i
deriv = diff(ydata)./diff(xdata);
    [~,ind]= min(deriv);
    center = xdata(ind);
    
    width = widthf(i);
    q = qf(i);
    
    H = Hf(i);
    A = 0;
    B = 0;
    C = 0;
%       w0    w    asymm
x0 = [center,width,q,H,A,B,C];
%for r/L=0.7 x0 = [center,width,q,0,0.3,2.25];
%Fano = a+bx+c*(q+2(x-x0)/w)^2/(1+4(x-x0)^2/w^2)
Fano = @(x,xdata)(x(7)+x(6)*xdata+x(5)*xdata.^2+x(4)*(x(3)+2*(xdata-x(1))/x(2)).^2./(1+4*(xdata-x(1)).^2/(x(2)^2)));

x=lsqcurvefit(Fano,x0,xdata,ydata);
  f=figure;
  plot(xdata,ydata,'-k',xdata,Fano(x,xdata),'-r');
  Q(i) = abs(x(1)/x(2));
end

N=length(data(:,1));
data = data(floor(2*N/5):N,:);
xdata=data(:,1)
Q_3 = zeros(10,1);
for i=1:10
    %i=2;

ydata=data(:,i+1)-min(data(:,i+1));
deriv = diff(ydata)./diff(xdata);
    [~,ind]= min(deriv);
    center = xdata(ind);
    
    width = widthf(i);
    q = qf(i);
    
    H = Hf(i);
    A = 0;
    B = 0;
    C = 0;
%       w0    w    asymm
x0 = [center,width,q,H,A,B,C];
%for r/L=0.7 x0 = [center,width,q,0,0.3,2.25];
%Fano = a+bx+c*(q+2(x-x0)/w)^2/(1+4(x-x0)^2/w^2)
Fano = @(x,xdata)(x(7)+x(6)*xdata+x(5)*xdata.^2+x(4)*(x(3)+2*(xdata-x(1))/x(2)).^2./(1+4*(xdata-x(1)).^2/(x(2)^2)));

x=lsqcurvefit(Fano,x0,xdata,ydata);
  f=figure;
  plot(xdata,ydata,'-k',xdata,Fano(x,xdata),'-r');
  Q_3(i) = abs(x(1)/x(2));
end
Q_3(1)=0;
theta_3=linspace(1.05,2,10);;

save('11_layers_1.05_to_3.0_new_results.mat','data','theta_3','Q_3')
plot(theta_3,Q_3,'-sg','LineWidth',2);
%end

