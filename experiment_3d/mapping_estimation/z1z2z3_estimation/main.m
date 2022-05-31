clear;
close all;

% load('premade_diff.mat')
% load('function/z1z2z3_estimation_func.mat')

tic


dk = 0.1;   % サンプル刻み
K_fin = 1.9;  % サンプリング終了時間
k = [0:dk:K_fin];

u1 = ones(length(k),1) * 0.5; % 並進速度
u2 = ones(length(k),1) * (0.6); % 回転角速度

s = zeros(length(k),3); % センサ変数 s = (s1,s2,s3) = (x,y,θ)
s(1,:) = [1 1 pi/4];    % 初期観測(初期位置)

s_corr = zeros(length(k),3); % 補正後のセンサ変数(z1,z3空間と等しい)、結果比較用
s_corr(1,:) = [0 0 0];


index_max = [2,11,2]; % 格子点インデックス{l,m,n}の最大値

index = zeros(length(k),3);  % 格子点インデックス{l,m,n}(>=1)

index_next = zeros(length(k),3);

index_real = zeros(length(k),3);  % 格子点のインデックス{l,m,n}の整数値

rho_tmp = zeros(length(k), index_max(1), index_max(2), index_max(3), 3);
rho = zeros(length(k),3);

m_case = zeros(length(k)-1, 1);  % 偏微分後関数選択のための場合分け


% 学習率
eta_s1 = 0; 
eta_s2 = 0;
% eta_s3 = 1.0 * 10 ^ (-3);

iteration = 1224; % 現時点で最適
% iteration = 2000;

E1_all_value = zeros(iteration-1,1);
Ereg_all_value = zeros(iteration-1,1);


grid_ = zeros(index_max(1), index_max(2), index_max(3), 3, iteration);

% 原点付近の格子点の固定
grid_(1,1,1,:,:) = [1 1 pi/4];
grid_(2,1,1,:,:) = [1+1/sqrt(2) 1+1/sqrt(2) pi/4];
grid_(1,2,1,:,:) = [1 1 pi/3];
grid_(1,1,2,:,:) = [1-1/sqrt(2) 1+1/sqrt(2) pi/4];


[s,s_corr] = sampling(s, s_corr, u1, u2, k, dk); % サンプル収集

grid_ = init_grid(grid_, index_max); % 格子点の初期値決定


% 最急降下法による格子点更新
for t = 1 : iteration

    break_switch = 0;
    
    grid_(:,:,:,:,t+1) = grid_(:,:,:,:,t);

    for j = 1 : length(k)

        % 格子点の選択
        [index, index_next, index_real, rho, rho_tmp, break_switch] = select_grid(grid_, t, j, index_max, index, index_next, index_real, rho, rho_tmp, break_switch);

        % 欠損時の補間
        if (break_switch == 0 && j > 1)

            index(j,:) = index(j-1,:);
            index_next(j,:) = index_next(j-1,:);     
            index_real(j,:) = index_real(j-1,:);
           
            [a,b,c] = index(j,:);

            rho(j,:) = rho_tmp(j,a,b,c,:);
            
            disp(j)
            disp("補間しました")

        end

        % 誤差関数の偏微分後関数選択のためのパターン分け
        if j > 1
            if index(j,2) == index(j-1,2)
                m_case(j-1) = 1;
            elseif index(j,2) == index_next(j-1,2) && index(j,2) ~= index(j-1,2)
                m_case(j-1) = 2;
            else
                m_case(j-1) = 3;
            end
        end

        break_switch = 0;

    end

    if t < iteration

        sd_update() % 最急降下法による格子点更新

        % 格子点の補助
        grid_(2,:,1,3,t+1) = grid_(1,:,1,3,t+1);
        grid_(1,:,2,3,t+1) = grid_(1,:,1,3,t+1);
        grid_(2,:,2,3,t+1) = grid_(1,:,1,3,t+1);

    else
        % z1,z2,z3の推定結果取得
        z = zeros(length(k),3);

        for j = 1:length(k)
            z(j,1) = index_real(j,1) + rho(j,1);
            z(j,2) = tan(pi/12) * (index_real(j,2) + rho(j,2));
            z(j,3) = index_real(j,3) + rho(j,3);
        end

    end

end    


%---g_t, h_tの生成--------------------------------------------------------

g_t = zeros(length(k)-1,1);
h_t = zeros(length(k)-1,1);

for j = 1 : length(k) - 1

    g_t(j) = ((z(j+1,2) - z(j,2)) * u1(j)) / ((z(j+1,1) - z(j,1)) * u2(j));
    h_t(j) = (z(j+1,1) - z(j,1)) / (u1(j) * dk);

end


disp("####################")



%--- z1, z2, z3の推定結果の描画 --------------------------------------

% 2D→1Dのグラフ（3D→1Dは可視化できないため）
figure;
hold on;
grid on;

axis([0.75 1.96 0.95 1.4 0 1.0]) % π/2 ≒ 1.57

plot3(s(:,3), s(:,1), s_corr(:,1), '--m', s(:,3), s(:,1), z1_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("s_3 = \theta [rad]",'fontsize',18)
ylabel("s_1 = x [m]",'fontsize',18)
zlabel("z_1",'fontsize',18)
legend(" True values: x'",' Estimated values','fontsize',20)

hold off;


figure;
hold on;
grid on;

axis([0.75 1.96 0.95 1.4 0 3.0]) % π/2 ≒ 1.57

plot3(s(:,3), s(:,1), tan(s_corr(:,3)), '--m', s(:,3), s(:,1), z2_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("s_3 = \theta [rad]",'fontsize',18)
ylabel("s_1 = x [m]",'fontsize',18)
zlabel("z_2",'fontsize',18)
legend(" True values: tan(\theta')",' Estimated values','fontsize',20)

hold off;


figure;
hold on;
grid on;

axis([0.75 1.96 0.95 1.4 0 1.0]) % π/2 ≒ 1.57

plot3(s(:,3), s(:,1), s_corr(:,2), '--m', s(:,3), s(:,1), z3_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("s_3 = \theta [rad]",'fontsize',18)
ylabel("s_1 = x [m]",'fontsize',18)
zlabel("z_3",'fontsize',18)
legend(" True values: y'",' Estimated values','fontsize',20)

hold off;


%---g_t, h_tの描画--------------------------------------------------------

figure;
hold on;
grid on;

axis([-0.2 1.4 -0.1 12]) % π/2 ≒ 1.57

g_ans = zeros(length(k)-1,1);

for j = 1 : length(k)-1
    g_ans(j) = 1 / (cos(s_corr(j,3)) * cos(s_corr(j,3)) * cos(s_corr(j,3)));
end

plot(s_corr(1:length(k)-1,3), g_ans(:), '--m', s_corr(1:length(k)-1,3), g_t(:),'o','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("\theta' [rad]",'fontsize',18)
ylabel('g_t', 'fontsize',18)
legend("真値：1/cos^3(\theta')",'推定値：g_t')

hold off;


figure;
hold on;
grid on;

axis([-0.2 1.4 0.2 1.4]) % π/2 ≒ 1.57

plot(s_corr(1:length(k)-1,3), cos(s_corr(1:length(k)-1,3)), '--m', s_corr(1:length(k)-1,3), h_t(:),'o','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("\theta' [rad]", 'fontsize',18)
ylabel("h_t", 'fontsize',18)
legend("真値：cos(\theta')",'推定値：h_t')

hold off;


%---格子点遷移の描画--------------------
figure;
hold on;
grid on;

axis([0 2 0.5 2.5 pi/6 pi])

for a = 1:2
    for b = 1:9
        for c = 1:2

            if b == 1
                plot3(grid_(a,1:9,c,1,1), grid_(a,1:9,c,2,1),grid_(a,1:9,c,3,1),'ko:')
            end
            if a == 1
                plot3(grid_(:,b,c,1,1), grid_(:,b,c,2,1),grid_(:,b,c,3,1),'k:')
            end
            
        end
    end
end

% plot3のための順番変更
for i = 1:2
    for a = 1:2
        for b = 1:9
            for d = 1:3
                tmp(i,a,b,d) = grid_(a,b,i,d,1);
            end
        end
    end
end

for a = 1:2
    for b = 1:9
        plot3(tmp(:,a,b,1),tmp(:,a,b,2),tmp(:,a,b,3),'k:');
    end
end

xlabel("x [m]",'fontsize',16)
ylabel("y [m]",'fontsize',16)
zlabel("θ [rad]",'fontsize',16)

hold off;


figure;
hold on;
grid on;

axis([0 2 0.5 2.5 pi/6 pi])

for a = 1:2
    for b = 1:9
        for c = 1:2

            if b == 1
                plot3(grid_(a,1:9,c,1,iteration), grid_(a,1:9,c,2,iteration),grid_(a,1:9,c,3,iteration),'ko:')
            end
            if a == 1
                plot3(grid_(:,b,c,1,iteration), grid_(:,b,c,2,iteration),grid_(:,b,c,3,iteration),'k:')
            end
            
        end
    end
end

% plot3のための順番変更
for i = 1:2
    for a = 1:2
        for b = 1:9
            for d = 1:3
                tmp2(i,a,b,d) = grid_(a,b,i,d,iteration);
            end
        end
    end
end

for a = 1:2
    for b = 1:9
        plot3(tmp2(:,a,b,1),tmp2(:,a,b,2),tmp2(:,a,b,3),'k:');
    end
end

xlabel("x [m]",'fontsize',16)
ylabel("y [m]",'fontsize',16)
zlabel("θ [rad]",'fontsize',16)

hold off;


%--- 誤差関数Eのグラフ --------------------------------------

figure;
hold on;
grid on;

axis([0 2000 -0.1 12]) % π/2 ≒ 1.57

plot(1:1:iteration-1, E1_all_value, '-k','MarkerEdgeColor','red','LineWidth', 1.5)
xlabel("t")
ylabel('E1')
% legend("真値：1/cos^3(s3')",'推定値：g')

hold off;


figure;
hold on;
grid on;

axis([0 2000 -0.01 0.1]) % π/2 ≒ 1.57

plot(1:1:iteration-1, Ereg_all_value, '-k','MarkerEdgeColor','red','LineWidth', 1.5)
xlabel("t")
ylabel('Ereg')
% legend("真値：1/cos^3(s3')",'推定値：g')

hold off;


% % matファイルへの保存
% save z1z2z3_estimation.mat

toc



% --- 以降関数 ----------------------------------------------

% function [s,s_corr] = sample(s, s_corr, u1, u2, k, dk)
%     for j = 1 : length(k) - 1
    
%         s(j+1,3) = s(j,3) + u2(j+1) * dk;
%         s(j+1,1) = s(j,1) + u1(j+1) * cos(s(j+1,3)) * dk;
%         s(j+1,2) = s(j,2) + u1(j+1) * sin(s(j+1,3)) * dk;
    
%         s_corr(j+1,3) = s_corr(j,3) + u2(j+1) * dk;
%         s_corr(j+1,1) = s_corr(j,1) + u1(j+1) * cos(s_corr(j+1,3)) * dk;
%         s_corr(j+1,2) = s_corr(j,2) + u1(j+1) * sin(s_corr(j+1,3)) * dk;
    
%     end
% end


f