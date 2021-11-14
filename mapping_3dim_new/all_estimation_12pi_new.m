clear;
close all;

tic

%---サンプル収集------------------------------------------------------------------------

dk1 = 0.1;   % 時間刻み
K1fin = 1.9;  %シミュレーション終了時間, length(k) = Kfin + 1
k1 = [0:dk1:K1fin];

u1_b1 = ones(length(k1),1) * 0.5; % 並進速度

if rem(i,1) == 1
    u2_b1 = ones(length(k1),1) * (0.6); % 回転角速度
else
    u2_b1 = ones(length(k1),1) * (0.6); % 回転角速度
end

si_b1 = zeros(length(k1),3); % 観測するセンサ変数 , s = (s1, s2, s3) = (x ,y, θ)
si_b1(1,:) = [1 1 pi/4];    % (s1, s2, s3)の初期値を設定

si_c1 = zeros(length(k1),3); % 補正後のセンサ変数(zi,z3空間と等しい)、結果比較用
si_c1(1,:) = [0 0 0];

for j = 1 : length(k1) - 1
    
    si_b1(j+1,3) = si_b1(j,3) + u2_b1(j+1) * dk1;
    si_b1(j+1,1) = si_b1(j,1) + u1_b1(j+1) * cos(si_b1(j+1,3)) * dk1;
    si_b1(j+1,2) = si_b1(j,2) + u1_b1(j+1) * sin(si_b1(j+1,3)) * dk1;

    si_c1(j+1,3) = si_c1(j,3) + u2_b1(j+1) * dk1;
    si_c1(j+1,1) = si_c1(j,1) + u1_b1(j+1) * cos(si_c1(j+1,3)) * dk1;
    si_c1(j+1,2) = si_c1(j,2) + u1_b1(j+1) * sin(si_c1(j+1,3)) * dk1;

end


%---原点付近の格子点探索・固定, およびその他初期値の決定(線形補間)------------------------------

l_max = 2;
m_max = 10;
n_max = 2;

iteration = 150;

param_s = zeros(l_max, m_max, n_max, 3, iteration);

param_s(1,1,1,:,1) = [1 1 pi/4];
param_s(2,1,1,:,1) = [1+1/sqrt(2) 1+1/sqrt(2) pi/4];
param_s(1,2,1,:,1) = [1 1 pi/3];
param_s(1,1,2,:,1) = [1-1/sqrt(2) 1+1/sqrt(2) pi/4];

s_l = param_s(2,1,1,:,1) - param_s(1,1,1,:,1);
s_m = param_s(1,2,1,:,1) - param_s(1,1,1,:,1);
s_n = param_s(1,1,2,:,1) - param_s(1,1,1,:,1);


for a = 1 : l_max
    for b = 1 : m_max
        for c = 1 : n_max

            if a <= l_max
                l_coef = a - 1;
            else
                l_coef = l_max - a;
            end

            if b <= m_max
                m_coef = b - 1;
            else
                m_coef = m_max - b;
            end

            if c <= n_max
                n_coef = c - 1;
            else
                n_coef = n_max - c;
            end

            param_s(a,b,c,:,1) = param_s(1,1,1,:,1) + l_coef * s_l + m_coef * s_m + n_coef * s_n;
    
        end
    end
end


%---偏微分後関数の生成------------------------------------------------------------------

%---E1の偏微分後関数の生成----------------

s = sym('s1',[4 4 4 3]); % l,m,n,iの順

for j = 1 : length(k1) - 1

    G = [s(2,1,1,:) - s(1,1,1,:); s(1,2,1,:) - s(1,1,1,:); s(1,1,2,:) - s(1,1,1,:)];
    H =  transpose(reshape(G,[3,3]));
    y = [si_b1(j,1) - s(1,1,1,1); si_b1(j,2) - s(1,1,1,2); si_b1(j,3) - s(1,1,1,3)];
    P = H \ y;

    for a = 1 : 3
        for b = 1 : 3
            for c = 1 : 3

                G2 = [s(a+1,b,c,:) - s(a,b,c,:); s(a,b+1,c,:) - s(a,b,c,:); s(a,b,c+1,:) - s(a,b,c,:)];
                H2 =  transpose(reshape(G2,[3,3]));
                y2 = [si_b1(j+1,1) - s(a,b,c,1); si_b1(j+1,2) - s(a,b,c,2); si_b1(j+1,3) - s(a,b,c,3)];
                P2 = H2 \ y2;

                e1 = P(2) - (P2(3) - P(3)) / (P2(1) - P(1));

                De1_s1_1{j,a,b,c} = matlabFunction(diff(e1,s(a,b,c,1)), 'vars', {s(:,:,:,:)});
                De1_s2_1{j,a,b,c} = matlabFunction(diff(e1,s(a,b,c,2)), 'vars', {s(:,:,:,:)});
                De1_s3_1{j,a,b,c} = matlabFunction(diff(e1,s(a,b,c,3)), 'vars', {s(:,:,:,:)});

                De1_s1_next_1{j,a+1,b,c} = matlabFunction(diff(e1,s(a,b,c,1)), 'vars', {s(:,:,:,:)});
                De1_s2_next_1{j,a+1,b,c} = matlabFunction(diff(e1,s(a,b,c,2)), 'vars', {s(:,:,:,:)});
                De1_s3_next_1{j,a+1,b,c} = matlabFunction(diff(e1,s(a,b,c,3)), 'vars', {s(:,:,:,:)});

                De1_s1_next_1{j,a,b+1,c} = matlabFunction(diff(e1,s(a,b,c,1)), 'vars', {s(:,:,:,:)});
                De1_s2_next_1{j,a,b+1,c} = matlabFunction(diff(e1,s(a,b,c,2)), 'vars', {s(:,:,:,:)});
                De1_s3_next_1{j,a,b+1,c} = matlabFunction(diff(e1,s(a,b,c,3)), 'vars', {s(:,:,:,:)});

                De1_s1_next_1{j,a,b,c+1} = matlabFunction(diff(e1,s(a,b,c,1)), 'vars', {s(:,:,:,:)});
                De1_s3_next_1{j,a,b,c+1} = matlabFunction(diff(e1,s(a,b,c,2)), 'vars', {s(:,:,:,:)});
                De1_s3_next_1{j,a,b,c+1} = matlabFunction(diff(e1,s(a,b,c,3)), 'vars', {s(:,:,:,:)});

            end
        end
    end

end

%---E4の偏微分後関数の生成----------------

E4 = 0;

for a = 1 : 2
    for b = 1 : m_max - 2
        for c = 1 : 2

            E4 = E4 + (( (s(a,b+2,c,1) - s(a,b+1,c,1)) ^ 2 + (s(a,b+2,c,2) - s(a,b+1,c,2)) ^ 2 + (s(a,b+2,c,3) - s(a,b+1,c,3)) ^ 2 )... 
                    - ( (s(a,b+1,c,1) - s(a,b,c,1)) ^ 2 + (s(a,b+1,c,2) - s(a,b,c,2)) ^ 2 + (s(a,b+1,c,3) - s(a,b,c,3)) ^ 2 )) ^ 2;
                
        end
    end
end



%---格子点選択および更新-----------------------------------------------------------------

imax = 10;

l_max_now = 2;
m_max_now = 4;
n_max_now = 2;

max_switch = 0;

for i = 1 : imax

    l_now = zeros(length(k1),1);
    m_now = zeros(length(k1),1);
    n_now = zeros(length(k1),1);

    l_next = zeros(length(k1),1);
    m_next = zeros(length(k1),1);
    n_next = zeros(length(k1),1);

    l_real = zeros(length(k1),1);
    m_real = zeros(length(k1),1);
    n_real = zeros(length(k1),1);

    rho_tmp = zeros(length(k1), l_max, m_max, n_max, 3);
    rho = zeros(length(k1),3);
    

    if i > 1

        param_s(:,:,:,:,1) = param_s(:,:,:,:,iteration);

        param_s(2,:,1,3,1) = param_s(1,:,1,3,1);
        param_s(1,:,2,3,1) = param_s(1,:,1,3,1);
        param_s(2,:,2,3,1) = param_s(1,:,1,3,1);

    end

    %---格子点の選択---------------------

    break_switch = 0;

    for j = 1 : length(k1)

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
        
                    A = [param_s(a2,b,c,:,1) - param_s(a,b,c,:,1); param_s(a,b2,c,:,1) - param_s(a,b,c,:,1); param_s(a,b,c2,:,1) - param_s(a,b,c,:,1)];

                    B = transpose(reshape(A,[3,3]));

                    x = [si_b1(j,1) - param_s(a,b,c,1,1); si_b1(j,2) - param_s(a,b,c,2,1); si_b1(j,3) - param_s(a,b,c,3,1)];

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

        % 誤差関数の偏微分後関数選択のための記憶
        if l_now(j+1) == l_now(j)
            a_mem(j) = 1;
        elseif l_now(j+1) == l_now(j) + 1
            a_mem(j) = 2;
        else
            a_mem(j) = 3;
        end

        if m_now(j+1) == m_now(j)
            b_mem(j) = 1;
        elseif m_now(j+1) == m_now(j) + 1
            b_mem(j) = 2;
        else
            b_mem(j) = 3;
        end

        if n_now(j+1) == n_now(j)
            c_mem(j) = 1;
        elseif n_now(j+1) == n_now(j) + 1
            c_mem(j) = 2;
        else
            c_mem(j) = 3;
        end


        break_switch = 0;

    end

    disp('-----')



    %---最急降下法による格子点更新-----------------------------------

    eta_s1 = 0.0 * 10 ^ (-6); % 学習率
    eta_s2 = 0.0 * 10 ^ (-6);
    eta_s3 = 1.0 * 10 ^ (-2);

    if i < 11
        iteration = 100;
    else
        iteration = 150;
    end

    if i < 4
        E4_coef = 70;
    else
        E4_coef = 35;
    end

    stop_switch = 0;

    % 格子点更新範囲の拡大
    if rem(i,1) == 0
        if l_max_now < l_max
            l_max_now = l_max_now + 1;
        else
            l_max_now = l_max;
        end

        if m_max_now < m_max
            m_max_now = m_max_now + 1;
        else
            m_max_now = m_max;
        end
    
        if n_max_now < n_max
            n_max_now = n_max_now + 1;
        else
            n_max_now = n_max;
        end
    end


    %---勾配法によるパラメータ更新---------------------------

    for t = 1 : iteration - 1

        param_s(:,:,:,:,t+1) = param_s(:,:,:,:,t);

        for j = 1 : length(k1)
            for a = 1 : l_max_now
                for b = 1 : m_max_now
                    for c = 1 : n_max_now

                        DE1_s1_2 = DE1_s1_1{a,b,c}(param_s(:,:,:,:,t));
                        DE1_s2_2 = DE1_s2_1{a,b,c}(param_s(:,:,:,:,t));
                        DE1_s3_2 = DE1_s3_1{a,b,c}(param_s(:,:,:,:,t));

                        param_s(a,b,c,1,t+1) = param_s(a,b,c,1,t) - eta_s1 * DE1_s1_2;
                        param_s(a,b,c,2,t+1) = param_s(a,b,c,2,t) - eta_s2 * DE1_s2_2;
                        param_s(a,b,c,3,t+1) = param_s(a,b,c,3,t) - eta_s3 * DE1_s3_2;

                        param_s(1,1,1,:,t+1) = [1 1 pi/4];
                        param_s(2,1,1,:,t+1) = [1+1/sqrt(2) 1+1/sqrt(2) pi/4];   
                        param_s(1,2,1,:,t+1) = [1 1 pi/3];
                        param_s(1,1,2,:,t+1) = [1-1/sqrt(2) 1+1/sqrt(2) pi/4];

                    end
                end
            end
        end

        % if t > 1

        %     if (E1_value(t) > E1_value(t-1))
        %         disp('t = ')
        %         disp(t)
        %         disp('Ef1が増加しました')
        %         disp('学習率を下げて再開します')
        %         disp('--------------------')
        %         eta_s1 = eta_s1 * 0.5; 
        %         eta_s2 = eta_s2 * 0.5;
        %         eta_s3 = eta_s3 * 0.5;
        %         param_s(:,:,:,:,t+1) = param_s(:,:,:,:,t);
        %         stop_switch = stop_switch + 1;
        %     end
        %     if stop_switch == 14
        %         disp('t = ')
        %         disp(t)
        %         disp('Ef1が増加しました')
        %         disp('iterationを強制終了します')
        %         disp("####################")
        %         iteration = t-1;
        %         break;
        %     end
        %     if t == iteration - 1
        %         iteration = t+1;
        %         disp('iterationを正常に終了することができました！')
        %         disp("####################")
        %         break
        %     end

        % end
    
    end


    if i > imax - 2

        param_s(2,:,1,3,iteration) = param_s(1,:,1,3,iteration);
        param_s(1,:,2,3,iteration) = param_s(1,:,1,3,iteration);
        param_s(2,:,2,3,iteration) = param_s(1,:,1,3,iteration);

        % z1,z2,z3の推定結果取得--------------------------------------

        z1_b1 = zeros(length(k1),1);
        z2_b1 = zeros(length(k1),1);
        z3_b1 = zeros(length(k1),1);

        l_now_2 = zeros(length(k1),1);
        m_now_2 = zeros(length(k1),1);
        n_now_2 = zeros(length(k1),1);

        l_next_2 = zeros(length(k1),1);
        m_next_2 = zeros(length(k1),1);
        n_next_2 = zeros(length(k1),1);

        l_real_2 = zeros(length(k1),1);
        m_real_2 = zeros(length(k1),1);
        n_real_2 = zeros(length(k1),1);

        rho_2_tmp = zeros(length(k1), l_max, m_max, n_max, 3);
        rho_2 = zeros(length(k1),3);


        for j = 1:length(k1)

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

        disp("####################")


        % 推定結果のplot--------------------------------------

        % figure;
        % hold on;
        % grid on;

        % axis([-5 5 -5 5]) % π/2 ≒ 1.57

        % plot(si_c1(:,3), si_c1(:,1), '--m', si_c1(:,3), z1_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z1 = f1(s) = s1 の答え合わせ
        % xlabel("s3' = θ")
        % ylabel('z1')
        % legend("真値：s1'",'推定値：z1')

        % hold off;

        figure;
        hold on;
        grid on;

        axis([-5 5 -5 5]) % π/2 ≒ 1.57

        plot(si_c1(:,3), tan(si_c1(:,3)), '--m', si_c1(:,3), z2_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z1 = f1(s) = s1 の答え合わせ
        xlabel("s3' = θ")
        ylabel("z2")
        legend("真値：tan(s3')",'推定値：z2')

        hold off;

        % figure;
        % hold on;
        % grid on;

        % axis([-5 5 -5 5]) % π/2 ≒ 1.57

        % plot(si_c1(:,3), si_c1(:,2), '--m', si_c1(:,3), z3_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z3 = f3(s) = s2 の答え合わせ
        % xlabel("s3' = θ")
        % ylabel('z3')
        % legend("真値：s2'",'推定値：z3')

        % hold off;



        %---g,hの導出--------------------------------------------------------

        g_b1 = zeros(length(k1),1);
        h_b1 = zeros(length(k1),1);

        for j = 1 : length(k1) - 1
        
            g_b1(j) = ((z2_b1(j+1) - z2_b1(j)) * u1_b1(j)) / ((z1_b1(j+1) - z1_b1(j)) * u2_b1(j));
            h_b1(j) = (z1_b1(j+1) - z1_b1(j)) / (u1_b1(j) * dk1);

        end

        g_b1(length(k1)) = 2 * g_b1(length(k1)-1) - g_b1(length(k1)-2);
        h_b1(length(k1)) = 2 * h_b1(length(k1)-1) - h_b1(length(k1)-2);


        % 推定結果のplot--------------------------------------

        figure;
        hold on;
        grid on;

        axis([-5 5 -5 5]) % π/2 ≒ 1.57

        g_ans = zeros(length(k1),1);

        for j = 1 : length(k1)
            g_ans(j) = 1 / (cos(si_c1(j,3)) * cos(si_c1(j,3)) * cos(si_c1(j,3)));
        end

        plot(si_c1(:,3), g_ans(:), '--m', si_c1(:,3), g_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
        xlabel("s3' = θ")
        ylabel('g')
        legend("真値：1/cos^3(s3')",'推定値：g')

        hold off;


        % figure;
        % hold on;
        % grid on;

        % axis([-5 5 -5 5]) % π/2 ≒ 1.57

        % plot(si_c1(:,3), cos(si_c1(:,3)), '--m', si_c1(:,3), h_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
        % xlabel("s3' = θ")
        % ylabel("h")
        % legend("真値：cos(s3')",'推定値：h')

        % hold off;

    end

end

toc