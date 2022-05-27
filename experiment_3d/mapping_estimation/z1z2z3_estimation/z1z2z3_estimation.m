clear;
close all;

load('function/z1z2z3_estimation_func.mat')

tic

% 時刻kで選択される格子点のインデックス{l,m,n}(>=1)
l_now = zeros(length(k),1);
m_now = zeros(length(k),1);
n_now = zeros(length(k),1);

% 時刻kで選択される格子点のインデックス{l,m,n} + 1
l_next = zeros(length(k),1);
m_next = zeros(length(k),1);
n_next = zeros(length(k),1);

% 時刻kで選択される格子点のインデックス{l,m,n}の整数値
l_real = zeros(length(k),1);
m_real = zeros(length(k),1);
n_real = zeros(length(k),1);

rho_tmp = zeros(length(k), l_max, m_max, n_max, 3);
rho = zeros(length(k),3);

% 偏微分後関数選択のためのパターン分け
b_mem = zeros(length(k)-1, 1); 


%---最急降下法による格子点更新-----------------------------------

% 学習率
eta_s1 = 0; 
eta_s2 = 0;
% eta_s3 = 1.0 * 10 ^ (-3);

iteration = 1224; % 現時点で最適
% iteration = 2000;

E1_all_value = zeros(iteration,1);
Ereg_all_value = zeros(iteration,1);


for t = 1 : iteration - 1

    %---格子点の選択---------------------
    break_switch = 0;
    
    param_s(:,:,:,:,t+1) = param_s(:,:,:,:,t);

    for j = 1 : length(k)

        for a = 1 : l_max - 1
            for b = 1 : m_max - 1
                for c = 1 : n_max - 1

                    if break_switch == 1 
                        break;
                    end
        
                    if a < l_max
                        a2 = a + 1;
                        l_real(j) = a - 1;
                    elseif a == l_max + 1
                        a2 = 1;
                        l_real(j) = -1;
                    else
                        a2 = a - 1;
                        l_real(j) = - a + l_max; 
                    end
        
                    if b < m_max
                        b2 = b + 1;
                        m_real(j) = b - 1;
                    elseif b == m_max + 1
                        b2 = 1;
                        m_real(j) = -1;
                    else
                        b2 = b - 1;
                        m_real(j) = - b + m_max;
                    end
        
                    if c < n_max
                        c2 = c + 1;
                        n_real(j) = c - 1;
                    elseif c == n_max + 1
                        c2 = 1;
                        n_real(j) = -1;
                    else
                        c2 = c - 1;
                        n_real(j) = - c + n_max;
                    end
        
                    A = [param_s(a2,b,c,:,t) - param_s(a,b,c,:,t); param_s(a,b2,c,:,t) - param_s(a,b,c,:,t); param_s(a,b,c2,:,t) - param_s(a,b,c,:,t)];

                    B = transpose(reshape(A,[3,3]));

                    x = [si_b1(j,1) - param_s(a,b,c,1,t); si_b1(j,2) - param_s(a,b,c,2,t); si_b1(j,3) - param_s(a,b,c,3,t)];

                    rho_tmp(j,a,b,c,:) = B \ x;

                    if (0 <= rho_tmp(j,a,b,c,1)) && (rho_tmp(j,a,b,c,1) <= 1)
                        if (0 <= rho_tmp(j,a,b,c,2)) && (rho_tmp(j,a,b,c,2) <= 1)
                            if (0 <= rho_tmp(j,a,b,c,3)) && (rho_tmp(j,a,b,c,3) <= 1)

                                rho(j,:) = rho_tmp(j,a,b,c,:);
        
                                l_now(j) = a;
                                m_now(j) = b;
                                n_now(j) = c;

                                l_next(j) = a2;
                                m_next(j) = b2;
                                n_next(j) = c2;

                                break_switch = 1;
        
                            end
                        end
                    end

                end
            end
        end

        % 欠損時の補間
        if (break_switch == 0 && j > 1)

            l_now(j) = l_now(j-1);
            m_now(j) = m_now(j-1);
            n_now(j) = n_now(j-1);

            l_next(j) = l_next(j-1);
            m_next(j) = m_next(j-1);
            n_next(j) = n_next(j-1);

            l_real(j) = l_real(j-1);
            m_real(j) = m_real(j-1);
            n_real(j) = n_real(j-1);

            a = l_now(j);
            b = m_now(j);
            c = n_now(j);

            rho(j,:) = rho_tmp(j,a,b,c,:);
            
            disp(j)
            disp("補間しました")

        end

        % 誤差関数の偏微分後関数選択のためのパターン分け
        if j > 1
            if m_now(j) == m_now(j-1)
                b_mem(j-1) = 1;
            elseif m_now(j) == m_next(j-1) && m_now(j) ~= m_now(j-1)
                b_mem(j-1) = 2;
            else
                b_mem(j-1) = 3;
            end
        end

        break_switch = 0;

    end

    DE1 = zeros(m_max,3);  % 誤差関数E_1の偏微分   
    DEreg = zeros(m_max,3);  % 正則化項の偏微分

    for b = 1 : m_max  % m = 1&2 はfix

        if b < 9
            eta_s3 = 7.5 * 10 ^ (-3);
            Ereg_coef = 1.0 * 10 ^ (-1);
        else
            eta_s3 = 7.5 * 10 ^ (-3);
            Ereg_coef = 2.0 * 10 ^ (-1);
        end

        for j = 1 : length(k) - 1
            
              % 時刻kの格子点とピッタリ一致した場合
              if b == m_now(j)

                if b_mem(j) == 1

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type1{j,1,1,1,x}(param_s(:,m_now(j):m_next(j),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type1_func{j}(param_s(:,m_now(j):m_next(j),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                
                elseif b_mem(j) == 2

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type2{j,1,1,1,x}(param_s(:,m_now(j):m_next(j+1),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type2_func{j}(param_s(:,m_now(j):m_next(j+1),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));

                else

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type3{j,1,1,1,x}(param_s(:,m_now(j):m_next(j),:,:,t), param_s(:,m_now(j+1):m_next(j+1),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type3_func{j}(param_s(:,m_now(j):m_next(j),:,:,t), param_s(:,m_now(j+1):m_next(j+1),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
            
                end

            elseif b == m_next(j) && b ~= m_now(j)
                
                if b_mem(j) == 1

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type1{j,1,2,1,x}(param_s(:,m_now(j):m_next(j),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type1_func{j}(param_s(:,m_now(j):m_next(j),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));

                elseif b_mem(j) == 2

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type2{j,1,2,1,x}(param_s(:,m_now(j):m_next(j+1),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type2_func{j}(param_s(:,m_now(j):m_next(j+1),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));

                else

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type3{j,1,2,1,x}(param_s(:,m_now(j):m_next(j),:,:,t), param_s(:,m_now(j+1):m_next(j+1),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type3_func{j}(param_s(:,m_now(j):m_next(j),:,:,t), param_s(:,m_now(j+1):m_next(j+1),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));

                end               
            
            elseif b == m_now(j+1) && b ~= m_now(j) && b ~= m_next(j)
                
                for x = 1 : 3
                    DE1(b,x) = DE1(b,x) + De1_type3{j,1,3,1,x}(param_s(:,m_now(j):m_next(j),:,:,t), param_s(:,m_now(j+1):m_next(j+1),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                end

                % 誤差関数E1の作成
                E1_all_value(t) = E1_all_value(t) + e1_type3_func{j}(param_s(:,m_now(j):m_next(j),:,:,t), param_s(:,m_now(j+1):m_next(j+1),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                
            elseif b == m_next(j+1) && b ~= m_now(j) && b ~= m_next(j) && b ~= m_now(j+1)
                
                if b_mem(j) == 2

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type2{j,1,3,1,x}(param_s(:,m_now(j):m_next(j+1),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type2_func{j}(param_s(:,m_now(j):m_next(j+1),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));

                else

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type3{j,1,4,1,x}(param_s(:,m_now(j):m_next(j),:,:,t), param_s(:,m_now(j+1):m_next(j+1),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type3_func{j}(param_s(:,m_now(j):m_next(j),:,:,t), param_s(:,m_now(j+1):m_next(j+1),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));

                end

            end
            
            if j == length(k) - 1

                % Eregの追加
                if b == 1                    
                    for x = 1 : 3
                        DEreg(b,x) = De_reg{1,1,1,x}(param_s(1,b:b+2,1,:,t));
                    end
                elseif b == m_max
                    for x = 1 : 3
                        DEreg(b,x) = De_reg{1,3,1,x}(param_s(1,b-2:b,1,:,t));
                    end
                elseif b == 2
                    for x = 1 : 3
                        DEreg(b,x) = De_reg{1,1,1,x}(param_s(1,b:b+2,1,:,t)) + De_reg{1,2,1,x}(param_s(1,b-1:b+1,1,:,t));
                    end
                elseif b == m_max - 1
                    for x = 1 : 3
                        DEreg(b,x) = De_reg{1,3,1,x}(param_s(1,b-2:b,1,:,t)) + De_reg{1,2,1,x}(param_s(1,b-1:b+1,1,:,t));
                    end
                else
                    for x = 1 : 3
                        DEreg(b,x) = De_reg{1,1,1,x}(param_s(1,b:b+2,1,:,t)) + De_reg{1,2,1,x}(param_s(1,b-1:b+1,1,:,t)) + De_reg{1,3,1,x}(param_s(1,b-2:b,1,:,t));
                    end
                end

                % 格子点の更新
                param_s(1,b,1,1,t+1) = param_s(1,b,1,1,t) - eta_s1 * DE1(b,1) - Ereg_coef * DEreg(b,1);
                param_s(1,b,1,2,t+1) = param_s(1,b,1,2,t) - eta_s2 * DE1(b,2) - Ereg_coef * DEreg(b,2);
                param_s(1,b,1,3,t+1) = param_s(1,b,1,3,t) - eta_s3 * DE1(b,3) - Ereg_coef * DEreg(b,3);

                param_s(1,1,1,:,t+1) = [1 1 pi/4];
                param_s(2,1,1,:,t+1) = [1+1/sqrt(2) 1+1/sqrt(2) pi/4];   
                param_s(1,2,1,:,t+1) = [1 1 pi/3];
                param_s(1,1,2,:,t+1) = [1-1/sqrt(2) 1+1/sqrt(2) pi/4];

            end

        end
        
    end

    param_s(2,:,1,3,t+1) = param_s(1,:,1,3,t+1);
    param_s(1,:,2,3,t+1) = param_s(1,:,1,3,t+1);
    param_s(2,:,2,3,t+1) = param_s(1,:,1,3,t+1);

    % 誤差関数Eregの作成
    for b = 2 : m_max-1
        Ereg_all_value(t+1) = Ereg_all_value(t+1) + e_reg_func(param_s(1,b-1:b+1,1,:,t+1));
    end

end


% 格子点の補助
% param_s(2,:,1,3,iteration) = param_s(1,:,1,3,iteration);
% param_s(1,:,2,3,iteration) = param_s(1,:,1,3,iteration);
% param_s(2,:,2,3,iteration) = param_s(1,:,1,3,iteration);

% z1,z2,z3の推定結果取得--------------------------------------

z1_b1 = zeros(length(k),1);
z2_b1 = zeros(length(k),1);
z3_b1 = zeros(length(k),1);

l_now_2 = zeros(length(k),1);
m_now_2 = zeros(length(k),1);
n_now_2 = zeros(length(k),1);

l_next_2 = zeros(length(k),1);
m_next_2 = zeros(length(k),1);
n_next_2 = zeros(length(k),1);

l_real_2 = zeros(length(k),1);
m_real_2 = zeros(length(k),1);
n_real_2 = zeros(length(k),1);

rho_2_tmp = zeros(length(k), l_max, m_max, n_max, 3);
rho_2 = zeros(length(k),3);


for j = 1:length(k)

    for a = 1 : l_max - 1
        for b = 1 : m_max - 1
            for c = 1 : n_max - 1

                if break_switch == 1 
                    break;
                end
    
                if a < l_max
                    a2 = a + 1;
                    l_real_2(j) = a - 1;
                elseif a == l_max + 1
                    a2 = 1;
                    l_real_2(j) = -1;
                else
                    a2 = a - 1;
                    l_real_2(j) = - a + l_max; 
                end
    
                if b < m_max
                    b2 = b + 1;
                    m_real_2(j) = b - 1;
                elseif b == m_max + 1
                    b2 = 1;
                    m_real_2(j) = -1;
                else
                    b2 = b - 1;
                    m_real_2(j) = - b + m_max;
                end
    
                if c < n_max
                    c2 = c + 1;
                    n_real_2(j) = c - 1;
                elseif c == n_max + 1
                    c2 = 1;
                    n_real_2(j) = -1;
                else
                    c2 = c - 1;
                    n_real_2(j) = - c + n_max;
                end
    
                A = [param_s(a2,b,c,:,iteration) - param_s(a,b,c,:,iteration); param_s(a,b2,c,:,iteration) - param_s(a,b,c,:,iteration); param_s(a,b,c2,:,iteration) - param_s(a,b,c,:,iteration)];

                B = transpose(reshape(A,[3,3]));

                x = [si_b1(j,1) - param_s(a,b,c,1,iteration); si_b1(j,2) - param_s(a,b,c,2,iteration); si_b1(j,3) - param_s(a,b,c,3,iteration)];

                rho_2_tmp(j,a,b,c,:) = B \ x;
                

                if (0 <= rho_2_tmp(j,a,b,c,1)) && (rho_2_tmp(j,a,b,c,1) <= 1)
                    if (0 <= rho_2_tmp(j,a,b,c,2)) && (rho_2_tmp(j,a,b,c,2) <= 1)
                        if (0 <= rho_2_tmp(j,a,b,c,3)) && (rho_2_tmp(j,a,b,c,3) <= 1)

                            rho_2(j,:) = rho_2_tmp(j,a,b,c,:);
    
                            l_now_2(j) = a;
                            m_now_2(j) = b;
                            n_now_2(j) = c;

                            l_next_2(j) = a2;
                            m_next_2(j) = b2;
                            n_next_2(j) = c2;

                            break_switch = 1;
    
                        end
                    end
                end

            end
        end
    end

    %欠損時の補間
    if (break_switch == 0) && (j > 1)

        l_now_2(j) = l_now_2(j-1);
        m_now_2(j) = m_now_2(j-1);
        n_now_2(j) = n_now_2(j-1);

        l_next_2(j) = l_next_2(j-1);
        m_next_2(j) = m_next_2(j-1);
        n_next_2(j) = n_next_2(j-1);

        l_real_2(j) = l_real_2(j-1);
        m_real_2(j) = m_real_2(j-1);
        n_real_2(j) = n_real_2(j-1);

        a = l_now_2(j);
        b = m_now_2(j);
        c = n_now_2(j);

        rho_2(j,:) = rho_2_tmp(j,a,b,c,:);
        
        disp(j)
        disp("補間しました")

    end

    break_switch = 0;

    z1_b1(j) = l_real_2(j) + rho_2(j,1);
    z2_b1(j) = tan(pi/12) * (m_real_2(j) + rho_2(j,2));
    z3_b1(j) = n_real_2(j) + rho_2(j,3);

end


%---g_t, h_tの生成--------------------------------------------------------

g_b1 = zeros(length(k)-1,1);
h_b1 = zeros(length(k)-1,1);

for j = 1 : length(k) - 1

    g_b1(j) = ((z2_b1(j+1) - z2_b1(j)) * u1_b1(j)) / ((z1_b1(j+1) - z1_b1(j)) * u2_b1(j));
    h_b1(j) = (z1_b1(j+1) - z1_b1(j)) / (u1_b1(j) * dk);

end


disp("####################")


%--- z1, z2, z3の推定結果の描画 --------------------------------------

% 2D→1Dのグラフ（3D→1Dは可視化できないため）
figure;
hold on;
grid on;

axis([0.75 1.96 0.95 1.4 0 1.0]) % π/2 ≒ 1.57

plot3(si_b1(:,3), si_b1(:,1), si_c1(:,1), '--m', si_b1(:,3), si_b1(:,1), z1_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("s_3 = \theta [rad]",'fontsize',18)
ylabel("s_1 = x [m]",'fontsize',18)
zlabel("z_1",'fontsize',18)
legend(" True values: x'",' Estimated values','fontsize',20)

hold off;


figure;
hold on;
grid on;

axis([0.75 1.96 0.95 1.4 0 3.0]) % π/2 ≒ 1.57

plot3(si_b1(:,3), si_b1(:,1), tan(si_c1(:,3)), '--m', si_b1(:,3), si_b1(:,1), z2_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("s_3 = \theta [rad]",'fontsize',18)
ylabel("s_1 = x [m]",'fontsize',18)
zlabel("z_2",'fontsize',18)
legend(" True values: tan(\theta')",' Estimated values','fontsize',20)

hold off;


figure;
hold on;
grid on;

axis([0.75 1.96 0.95 1.4 0 1.0]) % π/2 ≒ 1.57

plot3(si_b1(:,3), si_b1(:,1), si_c1(:,2), '--m', si_b1(:,3), si_b1(:,1), z3_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
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
    g_ans(j) = 1 / (cos(si_c1(j,3)) * cos(si_c1(j,3)) * cos(si_c1(j,3)));
end

plot(si_c1(1:length(k)-1,3), g_ans(:), '--m', si_c1(1:length(k)-1,3), g_b1(:),'o','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("\theta' [rad]",'fontsize',18)
ylabel('g_t', 'fontsize',18)
legend("真値：1/cos^3(\theta')",'推定値：g_t')

hold off;


figure;
hold on;
grid on;

axis([-0.2 1.4 0.2 1.4]) % π/2 ≒ 1.57

plot(si_c1(1:length(k)-1,3), cos(si_c1(1:length(k)-1,3)), '--m', si_c1(1:length(k)-1,3), h_b1(:),'o','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
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
                plot3(param_s(a,1:9,c,1,1), param_s(a,1:9,c,2,1),param_s(a,1:9,c,3,1),'ko:')
            end
            if a == 1
                plot3(param_s(:,b,c,1,1), param_s(:,b,c,2,1),param_s(:,b,c,3,1),'k:')
            end
            
        end
    end
end

% plot3のための順番変更
for i = 1:2
    for a = 1:2
        for b = 1:9
            for d = 1:3
                tmp(i,a,b,d) = param_s(a,b,i,d,1);
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
                plot3(param_s(a,1:9,c,1,iteration), param_s(a,1:9,c,2,iteration),param_s(a,1:9,c,3,iteration),'ko:')
            end
            if a == 1
                plot3(param_s(:,b,c,1,iteration), param_s(:,b,c,2,iteration),param_s(:,b,c,3,iteration),'k:')
            end
            
        end
    end
end

% plot3のための順番変更
for i = 1:2
    for a = 1:2
        for b = 1:9
            for d = 1:3
                tmp2(i,a,b,d) = param_s(a,b,i,d,iteration);
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

plot(0:1:iteration-1, E1_all_value, '-k','MarkerEdgeColor','red','LineWidth', 1.5)
xlabel("t")
ylabel('E1')
% legend("真値：1/cos^3(s3')",'推定値：g')

hold off;


figure;
hold on;
grid on;

axis([0 2000 -0.01 0.1]) % π/2 ≒ 1.57

plot(0:1:iteration-1, Ereg_all_value, '-k','MarkerEdgeColor','red','LineWidth', 1.5)
xlabel("t")
ylabel('Ereg')
% legend("真値：1/cos^3(s3')",'推定値：g')

hold off;


% matファイルへの保存
save z1z2z3_estimation.mat


toc