clear;
close all;

%--- サンプル収集-------------------------------------------
dk = 0.1;   % サンプル刻み
K_fin = 3.9;  % サンプリング終了時間
k = [0:dk:K_fin];

u1 = ones(length(k),1) * 0.25; % 並進速度
u2 = ones(length(k),1) * (0.3); % 回転角速度

s = zeros(length(k),3); % センサ変数 s = (s1,s2,s3) = (x,y,θ)
s(1,:) = [1 1 pi/4];    % 初期観測(初期位置)

s_corr = zeros(length(k),3); % 補正後のセンサ変数(z1,z3空間と等しい)、結果比較用
s_corr(1,:) = [0 0 0];

[s,s_corr] = sampling(s, s_corr, u1, u2, k, dk); % サンプル収集


hold on;
axis equal;
grid on;

axis([0 2 0.5 2.5])

xlabel("x",'FontSize',14)
ylabel("y",'FontSize',14)

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

plot(s(:,1),s(:,2),'-ko','MarkerEdgeColor','blue','MarkerFaceColor','blue', 'MarkerSize', 4)


