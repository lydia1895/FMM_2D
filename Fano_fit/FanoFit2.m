%function [ x ] = FanoFit2( data, width,q,H,A,B,C)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

load('13_layers_0.3_deg_right_peak.mat','data')
N=length(data(1,:))
xdata=data(:,1);%xdata=data(:,1);
for i=1:N-1
%i=2;
ydata=data(:,i+1)-min(data(:,i+1));
deriv = diff(ydata)./diff(xdata);
    [~,ind]= min(deriv);
    center = xdata(ind);
    
    width = 0.0003;
    q = 1.32;
    H = 0.03;
    A = 0;
    B = 0;%-13.3;
    C = 0;
%       w0    w    asymm
x0 = [center,width,q,H,A,B,C];
%for r/L=0.7 x0 = [center,width,q,0,0.3,2.25];
%Fano = a+bx+c*(q+2(x-x0)/w)^2/(1+4(x-x0)^2/w^2)
Fano = @(x,xdata)(x(7)+x(6)*xdata+x(5)*xdata.^2+x(4)*(x(3)+2*(xdata-x(1))/x(2)).^2./(1+4*(xdata-x(1)).^2/(x(2)^2)));

x=lsqcurvefit(Fano,x0,xdata,ydata);
  f=figure;
  plot(xdata,ydata,'-k',xdata,Fano(x,xdata),'-r');
  hold off
end
Q = x(1)/x(2)
save('13_layers_0.3_deg_right_peak.mat','data','Q')
%end

