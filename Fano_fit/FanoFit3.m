function [ x ] = FanoFit3( data, width_ar,q_ar,H_ar,A,B,C,min_val)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N=length(data(1,:))
xdata=data(:,1);%xdata=data(:,1);
for i=1:N-1
%i=2;
ydata=data(:,i+1)-min_val(i);%min(data(:,i+1));
deriv = diff(ydata)./diff(xdata);
    [~,ind]= min(deriv);
    center = xdata(ind);
%       w0    w    asymm
x0 = [center,width_ar(i),q_ar(i),H_ar(i),A,B,C];
%for r/L=0.7 x0 = [center,width,q,0,0.3,2.25];
%Fano = a+bx+c*(1+2(x-x0)/q/w)^2/(1+4(x-x0)^2/w^2)
Fano = @(x,xdata)(x(7)+x(6)*xdata+x(5)*xdata.^2+x(4)*(x(3)+2*(xdata-x(1))/x(2)).^2./(1+4*(xdata-x(1)).^2/(x(2)^2)));

x(i,:)=lsqcurvefit(Fano,x0,xdata,ydata);


  f=figure;
   plot(xdata,ydata,'-k',xdata,Fano(x(i,:),xdata),'or');

end
end



