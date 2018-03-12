
load('3_layers_025_to_075_results.mat','Q')
Q_1=Q;
load('3_layers_05_to_50_results.mat','Q')
Q_2 = Q;
theta_1 = [0.5 0.75]
theta_2 = linspace(1,5,9);
theta_full_3 = cat(2,theta_1,theta_2);
Q_full_3_layers=cat(2,Q_1(2:3),Q_2(2:10));
g3=figure
plot(theta_full_3,Q_full_3_layers,'-sr','Linewidth',2)

load('4_layers_05_to_50_results.mat','Q')
Q_1=Q;
load('4_layers_45_to_50_results.mat','Q')
Q_2 = Q;
theta_full_4 = linspace(0.5,5,10);
Q_full_4_layers=cat(2,Q_1(1:8),Q_2);
g4=figure
plot(theta_full_4,Q_full_4_layers,'-sg','Linewidth',2)

load('5_layers_05_to_15_results.mat','Q')
Q_1=Q;
load('5_layers_05_to_50_results.mat','Q')
Q_2 = Q;
theta_full_5 = linspace(0.5,5,10);
Q_full_5_layers=cat(2,Q_1(1:3),Q_2(4:10));
g5=figure
plot(theta_full_5,Q_full_5_layers,'-sb','Linewidth',2)


g6=figure
plot(theta_full_3,Q_full_3_layers,'-sr',theta_full_4,Q_full_4_layers,'-sg',...
    theta_full_5,Q_full_5_layers,'-sb','Linewidth',2)
h5 = legend('3 layers','4 layers','5 layers',3);
xlabel('theta, deg')
axis tight
ax = gca;
ax.XAxis.MinorTick = 'on';
ylabel('Q')
set(gca,'fontsize', 16)