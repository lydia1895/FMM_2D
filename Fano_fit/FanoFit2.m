%function [ x ] = FanoFit2( data, width,q,H,A,B,C)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%load('5_layers_05_to_15.mat','data')
%N=length(data(1,:))
N = length(R);

%thetamin = 0.5;
%thetamax = 5;
%theta_full = linspace(thetamin,thetamax,10);

xxdata=lambda;
left = 0.79;
right = 0.81;
[elem,num_elem]=min(abs(xxdata-left));
[elem2,num_elem_2]=min(abs(xxdata-right));



xdata = transpose(lambda(num_elem:num_elem_2));
ydata=R(num_elem:num_elem_2)-min(R(num_elem:num_elem_2));

deriv = diff(ydata)./diff(xdata);
    [~,ind]= min(deriv);
    center = xdata(ind);
    
    width = 0.001;
    q = -100;
    H = 0.1;
    A = 0;
    B = 0;%-13.3;
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
  hold off
  Q = abs(x(1)/x(2));


%g=figure
%plot(theta_full,Q,'-sg','Linewidth',2)
%theta(i) = theta_full(i);
%save('5_layers_05_to_15_results.mat','theta_full','Q','i')
%end
theta_Q = theta1(13);
%save('2_08_a_th_1.8_Q.mat','theta_Q','Q');

