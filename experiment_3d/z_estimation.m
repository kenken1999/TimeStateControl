clear;
close all;
tic

load('premade_diff.mat')


index_max = [2,9,2]; % 格子点インデックス{l,m,n}の最大値
index = zeros(length(k),3);  % 格子点インデックス{l,m,n}(>=1)
index_next = zeros(length(k),3);  % 格子点インデックス{l+1,m+1,n+1}
index_real = zeros(length(k),3);  % 格子点のインデックス{l,m,n}の実際の値

rho = zeros(length(k),3);

m_case = zeros(length(k)-1, 1);  % 偏微分後関数選択のための場合分け

iteration = 2000;

grid_ = zeros(index_max(1), index_max(2), index_max(3), 3, iteration);

% 原点付近の格子点の固定
for t = 1 : iteration
    grid_(1,1,1,:,t) = [1 1 pi/4];
    grid_(2,1,1,:,t) = [1+1/sqrt(2) 1+1/sqrt(2) pi/4];
    grid_(1,2,1,:,t) = [1 1 pi/4+atan(sigma(2))];
    grid_(1,1,2,:,t) = [1-1/sqrt(2) 1+1/sqrt(2) pi/4];
end

grid_ = init_grid(grid_, index_max); % 格子点の初期値決定

E1_all_value = zeros(iteration-1,1);
Ereg_all_value = zeros(iteration-1,1);
e_reg_value = zeros(iteration-1,1);
E_all_value = zeros(iteration-1,1);


%%%%% 最急降下法による格子点更新 %%%%%
for t = 1 : iteration
    grid_(:,:,:,:,t+1) = grid_(:,:,:,:,t);

    for j = 1 : length(k)
        % 格子点の選択
        [index, index_next, index_real, rho] = select_grid(grid_, s, t, j, index_max, index, index_next, index_real, rho);

        % 誤差関数の偏微分後関数選択のためのパターン分け
        if j > 1
            if index(j,2) == index(j-1,2)
                m_case(j-1) = 1;
            elseif (index(j,2) == index_next(j-1,2)) && (index(j,2) ~= index(j-1,2))
                m_case(j-1) = 2;
            else
                m_case(j-1) = 3;
            end
        end
    end

    if t < iteration
        % 最急降下法による格子点更新
        [grid_, E1_all_value, Ereg_all_value, E_all_value, e_reg_value] = sd_update(grid_, t, j, k, index_max, index, index_next, index_real, m_case, E1_all_value, Ereg_all_value, E_all_value, e_reg_value, De1_case1, De1_case2, De1_case3, e1_case1_func, e1_case2_func, e1_case3_func, e_reg_func, De_reg);

        % 格子点の補助
        grid_(2,:,1,3,t+1) = grid_(1,:,1,3,t+1);
        grid_(1,:,2,3,t+1) = grid_(1,:,1,3,t+1);
        grid_(2,:,2,3,t+1) = grid_(1,:,1,3,t+1);
    else
        % z1,z2,z3 の推定結果取得
        z = zeros(length(k),3);
        for j = 1:length(k)
            z(j,1) = sigma(1) * (index_real(j,1) + rho(j,1));
            z(j,2) = sigma(2) * (index_real(j,2) + rho(j,2));
            z(j,3) = sigma(3) * (index_real(j,3) + rho(j,3));
        end
    end
end    


%%%%% g_t, h_tの生成 %%%%%
g_t = zeros(length(k)-1,1);
h_t = zeros(length(k)-1,1);

for j = 1 : length(k) - 1
    g_t(j) = ((z(j+1,2) - z(j,2)) * u1(j)) / ((z(j+1,1) - z(j,1)) * u2(j));
    h_t(j) = (z(j+1,1) - z(j,1)) / (u1(j) * dk);
end


%%%%% z1, z2, z3の推定結果の描画 %%%%%
% figure;
% hold on;
% grid on;
% axis([pi/4-0.1 5*pi/8 0.95 1.4 0 1.0])
% plot3(s(:,3), s(:,1), s_corr(:,1), '--m', s(:,3), s(:,1), z(:,1),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
% xlabel("s_3 = \theta [rad]",'fontsize',18)
% ylabel("s_1 = x [m]",'fontsize',18)
% zlabel("z_1",'fontsize',18)
% legend(" True values: x'",' Estimated values','fontsize',20)
% hold off;

figure;
hold on;
grid on;
view([0,0])
axis([pi/4-0.1 5*pi/8 0.95 1.4 0 3.0])
plot3(s(:,3), s(:,1), tan(s_corr(:,3)), '--m', s(:,3), s(:,1), z(:,2),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("s_3 = \theta [rad]",'fontsize',18)
ylabel("s_1 = x [m]",'fontsize',18)
zlabel("z_2",'fontsize',18)
% legend("真値：tan(\theta')",'推定値')
% legend(" True values: tan(\theta')",' Estimated values','fontsize',20)
hold off;

% figure;
% hold on;
% grid on;
% axis([pi/4-0.1 5*pi/8 0.95 1.4 0 1.0])
% plot3(s(:,3), s(:,1), s_corr(:,2), '--m', s(:,3), s(:,1), z(:,3),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
% xlabel("s_3 = \theta [rad]",'fontsize',18)
% ylabel("s_1 = x [m]",'fontsize',18)
% zlabel("z_3",'fontsize',18)
% legend(" True values: y'",' Estimated values','fontsize',20)
% hold off;


%%%%% g_t, h_tの描画 %%%%%%
figure;
hold on;
grid on;
axis([pi/4-0.1 5*pi/8 -0.1 12])
g_def = zeros(length(k)-1,1);
for j = 1 : length(k)-1
    g_def(j) = 1 / (cos(s_corr(j,3)) * cos(s_corr(j,3)) * cos(s_corr(j,3)));
end
plot(s(1:length(k)-1,3), g_def(:), '--m', s(1:length(k)-1,3), g_t(:),'o','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("s_3 = \theta [rad]",'fontsize',18)
ylabel('g_t', 'fontsize',18)
legend("真値：1/cos^3(\theta')",'推定値：g_t')
hold off;

% figure;
% hold on;
% grid on;
% axis([pi/4-0.1 5*pi/8 0.2 1.4])
% plot(s(1:length(k)-1,3), cos(s_corr(1:length(k)-1,3)), '--m', s(1:length(k)-1,3), h_t(:),'o','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
% xlabel("s_3 = \theta [rad]", 'fontsize',18)
% ylabel("h_t", 'fontsize',18)
% legend("真値：cos(\theta')",'推定値：h_t')
% hold off;


%%%%% 格子点遷移の描画 %%%%%%
% figure;
% hold on;
% grid on;
% axis([0 2 0.5 2.5 pi/6 pi])
% for l = 1:index_max(1)
%     for m = 1:index_max(2)
%         for n = 1:index_max(3)
%             if m == 1
%                 plot3(grid_(l,1:index_max(2),n,1,1), grid_(l,1:index_max(2),n,2,1),grid_(l,1:index_max(2),n,3,1),'ko:')
%             end
%             if l == 1
%                 plot3(grid_(:,m,n,1,1), grid_(:,m,n,2,1),grid_(:,m,n,3,1),'k:')
%             end         
%         end
%     end
% end

% % plot3のための順番変更
% for i = 1:2
%     for l = 1:2
%         for m = 1:index_max(2)
%             for d = 1:3
%                 tmp(i,l,m,d) = grid_(l,m,i,d,1);
%             end
%         end
%     end
% end

% for l = 1:2
%     for m = 1:index_max(2)
%         plot3(tmp(:,l,m,1),tmp(:,l,m,2),tmp(:,l,m,3),'k:');
%     end
% end

% plot3(s(:,1),s(:,2),s(:,3),'-bo','MarkerEdgeColor','blue','MarkerFaceColor','blue','MarkerSize', 4);
% xlabel("x [m]",'fontsize',16)
% ylabel("y [m]",'fontsize',16)
% zlabel("θ [rad]",'fontsize',16)
% hold off;


figure;
hold on;
grid on;
view([0,0])
axis([0 2 0.5 2.5 pi/6 pi])
for l = 1:index_max(1)
    for m = 1:index_max(2)
        for n = 1:index_max(3)
            if m == 1
                plot3(grid_(l,1:index_max(2),n,1,iteration), grid_(l,1:index_max(2),n,2,iteration),grid_(l,1:index_max(2),n,3,iteration),'ko:')
            end
            if l == 1
                plot3(grid_(:,m,n,1,iteration), grid_(:,m,n,2,iteration),grid_(:,m,n,3,iteration),'k:')
            end           
        end
    end
end

% plot3のための順番変更
for i = 1:2
    for l = 1:2
        for m = 1:index_max(2)
            for d = 1:3
                tmp2(i,l,m,d) = grid_(l,m,i,d,iteration);
            end
        end
    end
end

for l = 1:2
    for m = 1:index_max(2)
        plot3(tmp2(:,l,m,1),tmp2(:,l,m,2),tmp2(:,l,m,3),'k:');
    end
end

plot3(s(:,1),s(:,2),s(:,3),'-bo','MarkerEdgeColor','blue','MarkerFaceColor','blue','MarkerSize', 4);
xlabel("x [m]",'fontsize',16)
ylabel("y [m]",'fontsize',16)
zlabel("θ [rad]",'fontsize',16)
hold off;


%%%%% 誤差関数 E_1, E_reg, e_reg, E の値 %%%%%%
figure;
hold on;
grid on;
axis([0 iteration -0.1 12])
plot(0:1:iteration-2, E1_all_value, '-k','MarkerEdgeColor','red','LineWidth', 1.5)
xlabel("t")
ylabel('E1')
hold off;

% figure;
% hold on;
% grid on;
% axis([0 iteration -0.01 0.1])
% plot(0:1:iteration-2, Ereg_all_value, '-k','MarkerEdgeColor','red','LineWidth', 1.5)
% xlabel("t")
% ylabel('Ereg')
% hold off;
% figure;
% hold on;
% grid on;

% axis([0 iteration -1.0*10^(-5) 1.0*10^(-3)]) 
% plot(0:1:iteration-2, e_reg_value, '-k','MarkerEdgeColor','red','LineWidth', 1.5)
% xlabel("t")
% ylabel('e_reg')
% hold off;

% figure;
% hold on;
% grid on;
% axis([0 iteration -0.01 7])
% plot(0:1:iteration-2, E_all_value, '-k','MarkerEdgeColor','red','LineWidth', 1.5)
% xlabel("t")
% ylabel('E')
% hold off;


%%%%% matファイルへの保存 %%%%%
save z_estimation.mat


toc