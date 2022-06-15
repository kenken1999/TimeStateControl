clear;
close all;
load('../z1z2z3_estimation/main.mat')

tic


imax_g = 5;
imax_h = 15;

mean_g = zeros(imax_g,3);
var_g = 0.5;

mean_h = zeros(imax_h,3);
var_h = 0.5;

for i = 1 : imax_g
    mean_g(i,1) = i * 0.5;
    mean_g(i,2) = i * 0.5;
    mean_g(i,3) = i * 0.5;
end

for i = 1 : imax_h
    mean_h(i,1) = i * 0.5;
    mean_h(i,2) = i * 0.5;
    mean_h(i,3) = i * 0.5;
end


phi_g = zeros(length(k)-1,imax_g);
phi_h = zeros(length(k)-1,imax_h);

w_g = zeros(imax_g,1);
w_h = zeros(imax_h,1);

g_est = zeros(length(k)-1,1);
h_est = zeros(length(k)-1,1);


for j = 1 : length(k)-1

    for i = 1 : imax_g
        phi_g(j,i) = exp( - sqrt( (s(j,1) - mean_g(i,1)) ^ 2 + (s(j,2) - mean_g(i,2)) ^ 2 + (s(j,3) - mean_g(i,3)) ^ 2 ) / (2 * var_g ^ 2) );
    end
    
    for i = 1 : imax_h
        phi_h(j,i) = exp( - sqrt( (s(j,1) - mean_h(i,1)) ^ 2 + (s(j,2) - mean_h(i,2)) ^ 2 + (s(j,3) - mean_h(i,3)) ^ 2 ) / (2 * var_h ^ 2) );
    end

end


w_g = (transpose(phi_g) * phi_g) \ transpose(phi_g) * g_t;
w_h = (transpose(phi_h) * phi_h) \ transpose(phi_h) * h_t;


g_est = phi_g * w_g;
h_est = phi_h * w_h;


%---g,hのplot--------------------------------------------------------

figure;
hold on;
grid on;

axis([0.75 1.96 0.95 1.4 0 15.0]) % π/2 ≒ 1.57

for i = 1 : length(k)-1
    g_ans(i) = 1 / (cos(s_corr(i,3)) * cos(s_corr(i,3)) * cos(s_corr(i,3)));
end

plot3(s(1:length(k)-1,3),s(1:length(k)-1,1), g_ans(:), '--m', s(1:length(k)-1,3),s(1:length(k)-1,1), g_est(:), '-k', s(1:length(k)-1,3), s(1:length(k)-1,1), g_t(:),'o','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("s_3 = \theta [rad]",'fontsize',18)
ylabel("s_1 = x [m]",'fontsize',18)
zlabel('g','FontSize',18)
legend(" True mapping： 1/cos^3(\theta')",' Estimated mapping： g',' Training data','FontSize',18)

hold off;


figure;
hold on;
grid on;

axis([0.75 1.96 0.95 1.4 0.4 1.2]) % π/2 ≒ 1.57

plot3(s(1:length(k)-1,3),s(1:length(k)-1,1), cos(s_corr(1:length(k)-1,3)), '--m', s(1:length(k)-1,3), s(1:length(k)-1,1),h_est(:), '-k', s(1:length(k)-1,3),s(1:length(k)-1,1), h_t(:),'o','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("s_3 = \theta [rad]",'fontsize',18)
ylabel("s_1 = x [m]",'fontsize',18)
zlabel("h",'FontSize',18)
legend(" True mapping： cos(\theta')",' Estimated mapping： h',' Training data','FontSize',18)

hold off;


% matファイルへの保存
save gh_estimation.mat

toc