clear;
close all;
load('all_estimation_12pi_new.mat')

tic

%---サンプル収集と誤差関数の定義----------------------------------------------------

dk1 = 0.1;   % 時間刻み
K1fin = 1.9;  %シミュレーション終了時間, length(k) = Kfin + 1
k1 = [0:dk1:K1fin];

u1_b1 = ones(length(k1),1) * 0.5; % 並進速度
u2_b1 = ones(length(k1),1) * (0.6); % 回転角速度

si_b1 = zeros(length(k1),3); % 観測するセンサ変数 , s = (s1, s2, s3) = (x ,y, θ)
si_b1(1,:) = [1 1 pi/4];    % (s1, s2, s3)の初期値を設定

si_c1 = zeros(length(k1),3); % 補正後のセンサ変数(zi,z3空間と等しい)、結果比較用
si_c1(1,:) = [0 0 0];


imax_g = 20;
imax_h = 10;

mean_g = zeros(imax_g,3);
var_g = 1;

mean_h = zeros(imax_h,3);
var_h = 1;

for i = 1 : imax_g
    mean_g(i,1) = i * 0.1;
    mean_g(i,2) = i * 0.1;
    mean_g(i,3) = i * 0.1;
end

for i = 1 : imax_h
    mean_h(i,1) = i * 0.1 - 0.25;
    mean_h(i,2) = i * 0.1 - 0.25;
    mean_h(i,3) = i * 0.1 - 0.25;
end


phi_g = zeros(length(k1)-1,imax_g);
phi_h = zeros(length(k1)-1,imax_h);

w_g = zeros(imax_g,1);
w_h = zeros(imax_h,1);

g_est = zeros(length(k1)-1,1);
h_est = zeros(length(k1)-1,1);


for j = 1 : length(k1)-1

    for i = 1 : imax_g
        phi_g(j,i) = exp( - sqrt( (si_b1(j,1) - mean_g(i,1)) ^ 2 + (si_b1(j,2) - mean_g(i,2)) ^ 2 + (si_b1(j,3) - mean_g(i,3)) ^ 2 ) / (2 * var_g ^ 2) );
    end
    
    for i = 1 : imax_h
        phi_h(j,i) = exp( - sqrt( (si_b1(j,1) - mean_h(i,1)) ^ 2 + (si_b1(j,2) - mean_h(i,2)) ^ 2 + (si_b1(j,3) - mean_h(i,3)) ^ 2 ) / (2 * var_h ^ 2) );
    end

    if j < length(k1) - 1

        si_b1(j+1,3) = si_b1(j,3) + u2_b1(j+1) * dk1;
        si_b1(j+1,1) = si_b1(j,1) + u1_b1(j+1) * cos(si_b1(j+1,3)) * dk1;
        si_b1(j+1,2) = si_b1(j,2) + u1_b1(j+1) * sin(si_b1(j+1,3)) * dk1;

        si_c1(j+1,3) = si_c1(j,3) + u2_b1(j+1) * dk1;
        si_c1(j+1,1) = si_c1(j,1) + u1_b1(j+1) * cos(si_c1(j+1,3)) * dk1;
        si_c1(j+1,2) = si_c1(j,2) + u1_b1(j+1) * sin(si_c1(j+1,3)) * dk1;

    end

end


w_g = (transpose(phi_g) * phi_g) \ transpose(phi_g) * g_b1;
w_h = (transpose(phi_h) * phi_h) \ transpose(phi_h) * h_b1;


g_est = phi_g * w_g;
h_est = phi_h * w_h;


%---g,hのplot--------------------------------------------------------

figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

for i = 1 : length(k1)-1
    g_ans(i) = 1 / (cos(si_c1(i,3)) * cos(si_c1(i,3)) * cos(si_c1(i,3)));
end

plot(si_c1(1:length(k1)-1,3), g_ans(:), '--m', si_c1(1:length(k1)-1,3), g_est(:), '-k', si_c1(1:length(k1)-1,3), g_b1(:),'o','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("s3' = θ")
ylabel('g')
legend("真値：1/cos^3(s3')",'推定値：g')

hold off;


figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c1(1:length(k1)-1,3), cos(si_c1(1:length(k1)-1,3)), '--m', si_c1(1:length(k1)-1,3), h_est(:), '-k', si_c1(1:length(k1)-1,3), h_b1(:),'o','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("s3' = θ")
ylabel("h")
legend("真値：cos(s3')",'推定値：h')

hold off;



toc