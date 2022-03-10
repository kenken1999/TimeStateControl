clear;
close all;
load('all_estimation_12pi_new.mat')


hold on;
axis equal;
grid on;

axis([0 2 0.5 2.5])

xlabel("x",'FontSize',14)
ylabel("y",'FontSize',14)


% plot(1+1/4,2,'kx','MarkerSize', 15,'LineWidth',4)
% plot(1,1,'rx','MarkerSize', 15,'LineWidth',4)

z1_plot = 0:1:3;
z3_plot = 0:1:3;

for j = 0 : 2
    if j == 1
        plot(z1_plot, z1_plot - sqrt(2) + sqrt(2) * j, '-r')
        plot(z3_plot, - z3_plot + 2 + sqrt(2) * (j-1), '-b')
    else
        plot(z1_plot, z1_plot - sqrt(2) + sqrt(2) * j, '--k')
        plot(z3_plot, - z3_plot + 2 + sqrt(2) * (j-1), '--k')
    end
end

plot(si_b1(:,1),si_b1(:,2),'-ko','MarkerEdgeColor','blue','MarkerFaceColor','blue', 'MarkerSize', 4)


