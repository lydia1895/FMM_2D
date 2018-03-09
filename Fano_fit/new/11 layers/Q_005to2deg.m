


load('11_layers_0.05_to_0.50_deg.mat','Q')
load('11_layers_0.55_to_1.0_new.mat','Q_2')
load('11_layers_1.05_to_3.0_new_results.mat','Q_3')

Q=cat(1,Q,Q_2,Q_3);
for i=1:30
    if Q(i)>10000
        Q(i)=0
    end
end
theta = linspace(0.05,1,20)
theta_2 = linspace(1.05,2,10)
theta = cat(2,theta,theta_2);
f=figure;
    plot(theta,Q, '-sg', 'LineWidth', 2);
    axis tight
    ax = gca;
    ax.XAxis.MinorTick = 'on';
    xlabel('theta,deg')
    ylabel('Q')
    set(gca,'fontsize', 16)