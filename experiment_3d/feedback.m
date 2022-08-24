clear;
close all;
tic

load('z_estimation.mat')
load('gh_estimation.mat')


dt = 0.01;  % 観測間隔
T_fin = 3;  % シミュレーション終了時間
t = [0:dt:T_fin];

s_f = zeros(length(t),3);
s_f(1,:) = [1.25 2 pi/2];  % 初期観測（初期位置）

z_f = zeros(length(t),3);

u1_f = ones(length(t),1) * (-1);
u2_f = ones(length(t),1) * (1);

mu1 = ones(length(t),1) * (-1);  % μ1 = v1 = u1cosθ
mu2 = ones(length(t),1);  % μ2 = v2/v1 = (u2/u1)*(1/(cosθ)^3)

g_map = zeros(1,length(t));
h_map = zeros(1,length(t));

index_f = zeros(length(t),3);
index_f_next = zeros(length(t),3);
index_f_real = zeros(length(t),3);
rho_f = zeros(length(t),3);

phi_g_now = zeros(1,number_g);
phi_h_now = zeros(1,number_h);

% feedback係数
k2 = 4;
k3 = 5;


%%%%% 制御範囲の描画 %%%%%
hold on;
axis equal;
grid on;
axis([0 2 0.5 2.5])
xlabel("x",'FontSize',14)
ylabel("y",'FontSize',14)

% 初期位置と原点の描画
plot(s_f(1,1), s_f(1,2), 'kx', 'MarkerSize', 10, 'LineWidth', 2)
plot(1, 1, 'rx', 'MarkerSize', 10, 'LineWidth', 2)

% z1-z3座標の描画
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


%%%%% 学習結果を用いたFeedback制御 %%%%%
for i = 1:length(t)-1
    %%% 制御時のアニメーションを表示する場合、以下を使用 %%%
    if i == 1
        h = plot(s_f(1,1), s_f(1,3), 'o', 'MarkerSize' ,20, 'MarkerFaceColor', 'b');
        h2 = plot(s_f(1,1) + 0.1*cos(s_f(1,3)), s_f(1,2) + 0.1*sin(s_f(1,3)), 'o', 'MarkerSize' ,8, 'MarkerFaceColor', 'r');
    end
    set(h, 'XData', s_f(i,1),'YData', s_f(i,2));
    set(h2, 'XData', s_f(i,1) + 0.1*cos(s_f(i,3)),'YData', s_f(i,2) + 0.1*sin(s_f(i,3)));
    set(gcf, 'Color', 'white');  % figureの背景を透明に設定
    set(gca, 'Color', 'white');  % axisの背景を透明に設定
    drawnow;
    xlabel("x", 'FontSize', 16)
    ylabel("y", 'FontSize', 16)

    %%% gifを保存する場合、以下を使用 %%%
    F = getframe(gcf);
    [X,map] = rgb2ind(F.cdata, 256);  % RGBデータをインデックス付きデータに変更
    if i == 1
        imwrite(X,map, 'feedback.gif')  % GIFファイルに書き出し
    else
        imwrite(X,map, 'feedback.gif', 'WriteMode', 'append')  % 2回目以降は'append'でアニメーションを作成
    end

    %%% 車両の軌跡を表示する場合、以下を使用 %%%
    % hold on;
    % plot(s_f(i,1),s_f(i,2), 'o', 'MarkerSize' ,3, 'MarkerFaceColor', 'b','MarkerEdgeColor', 'b');

    % 格子点の選択
    [index_f, index_f_next, index_f_real, rho_f] = select_grid(grid_, s_f, iteration, i, index_max, index_f, index_f_next, index_f_real, rho_f);

    % z1,z2,z3 の出力
    z_f(i,1) = sigma(1) * (index_f_real(i,1) + rho_f(i,1));
    z_f(i,2) = sigma(2) * (index_f_real(i,2) + rho_f(i,2));
    z_f(i,3) = sigma(3) * (index_f_real(i,3) + rho_f(i,3));

    for p = 1 : number_g
        phi_g_now(p) = exp( - ((s_f(i,1) - mean_g(p,1))^2 + (s_f(i,2) - mean_g(p,2))^2 + (s_f(i,3) - mean_g(p,3))^2) / (2 * var_g^2));
    end
    
    for q = 1 : number_h
        phi_h_now(q) = exp( - ((s_f(i,1) - mean_h(q,1))^2 + (s_f(i,2) - mean_h(q,2))^2 + (s_f(i,3) - mean_h(q,3))^2) / (2 * var_h^2));
    end

    g_map(i) = phi_g_now * w_g;
    h_map(i) = phi_h_now * w_h;

    mu1(i) = -2 * z_f(i,1);  % μ1 = -λz1 (λは正の定数)
     
    % μ1の正負に応じてμ2を場合分け
    if mu1(i) > 0
        mu2(i) = -k2 * z_f(i,2) - k3 * z_f(i,3);
    else
        mu2(i) = k2 * z_f(i,2) - k3 * z_f(i,3); 
    end

    u1_f(i) = mu1(i) / h_map(i);
    u2_f(i) = mu2(i) * u1_f(i) / g_map(i);

    % 次時刻の観測
    s_f(i+1,3) = s_f(i,3) + u2_f(i) * dt;
    s_f(i+1,1) = s_f(i,1) + u1_f(i) * cos(s_f(i+1,3)) * dt;
    s_f(i+1,2) = s_f(i,2) + u1_f(i) * sin(s_f(i+1,3)) * dt;
end


%%%%% (x,y,θ) の Feedback結果の図を表示する場合、以下を使用 %%%%%
% target = ones(length(t),1);
% target2 = ones(length(t),1) * (pi/4);
% tiledlayout(3,1);
% nexttile
% plot(t, target, '--k', t, s_f(:,1), '-k','LineWidth', 2)
% ylabel("x [m]",'fontsize',16)
% nexttile
% plot(t, target, '--k', t, s_f(:,2), '-k','LineWidth', 2)
% ylabel("y [m]",'fontsize',16)
% nexttile
% plot(t, target2, '--k', t, s_f(:,3), '-k','LineWidth', 2)
% xlabel("time [s]",'fontsize',16)
% ylabel("θ [rad]",'fontsize',16)


toc