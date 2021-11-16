clear;
close all;

load('all_estimation_12pi_func.mat')

tic


%---格子点選択および更新-----------------------------------------------------------------

imax = 20;

l_max_now = 2;
m_max_now = 4;
n_max_now = 2;

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

    b_mem = zeros(length(k1)-1, 1);
    

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

        % 誤差関数の偏微分後関数選択のための分類
        if j > 1
            if m_now(j) == m_now(j-1)
                b_mem(j-1) = 1;
            elseif m_now(j) == m_next(j-1)
                b_mem(j-1) = 2;
            else
                b_mem(j-1) = 3;
            end
        end

        break_switch = 0;

    end

    disp('-----')


    %---最急降下法による格子点更新-----------------------------------

    eta_s1 = 0.0 * 10 ^ (-6); % 学習率
    eta_s2 = 0.0 * 10 ^ (-6);
    eta_s3 = 1.0 * 10 ^ (-4);

    if i < 11
        iteration = 100;
    elseif i > 18
        iteration = 150;
    else
        iteration = 100;
    end

    if i < 4
        Ereg_coef = 50;
    else
        Ereg_coef = 50;
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


    for t = 1 : iteration - 1

        param_s(:,:,:,:,t+1) = param_s(:,:,:,:,t);

        DE1 = zeros(m_max,3);
        DEreg = zeros(m_max,3);

        for b = 3 : m_max_now % m = 1&2 はfix

            for j = 1 : length(k1) - 1 
              
                % 時刻kの格子点とピッタリ一致した場合
                if b == m_now(j)

                    if b_mem(j) == 1
                        for x = 1 : 3
                            DE1(b,x) = DE1(b,x) + De1_type1{j,1,1,1,x}(param_s(:,m_now(j):m_next(j),:,:,t));
                        end
                    elseif b_mem(j) == 2
                        for x = 1 : 3
                            DE1(b,x) = DE1(b,x) + De1_type2{j,1,1,1,x}(param_s(:,m_now(j):m_next(j+1),:,:,t));
                        end
                    else
                        for x = 1 : 3
                            DE1(b,x) = DE1(b,x) + De1_type3{j,1,1,1,x}(param_s(:,m_now(j):m_next(j),:,:,t), param_s(:,m_now(j+1):m_next(j+1),:,:,t));
                        end
                    end

                elseif b == m_next(j) && b ~= m_now(j)
                    
                    if b_mem(j) == 1
                        for x = 1 : 3
                            DE1(b,x) = DE1(b,x) + De1_type1{j,1,2,1,x}(param_s(:,m_now(j):m_next(j),:,:,t));
                        end
                    elseif b_mem(j) == 2
                        for x = 1 : 3
                            DE1(b,x) = DE1(b,x) + De1_type2{j,1,2,1,x}(param_s(:,m_now(j):m_next(j+1),:,:,t));
                        end
                    else
                        for x = 1 : 3
                            DE1(b,x) = DE1(b,x) + De1_type3{j,1,2,1,x}(param_s(:,m_now(j):m_next(j),:,:,t), param_s(:,m_now(j+1):m_next(j+1),:,:,t));
                        end
                    end               
                
                elseif b == m_now(j+1) && b ~= m_now(j) && b ~= m_next(j)
                    
                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type3{j,1,3,1,x}(param_s(:,m_now(j):m_next(j),:,:,t), param_s(:,m_now(j+1):m_next(j+1),:,:,t));
                    end
                    
                elseif b == m_next(j+1) && b ~= m_now(j) && b ~= m_next(j) && b == m_now(j+1)
                    
                    if b_mem(j) == 2
                        for x = 1 : 3
                            DE1(b,x) = DE1(b,x) + De1_type2{j,1,3,1,x}(param_s(:,m_now(j):m_next(j+1),:,:,t));
                        end
                    else
                        for x = 1 : 3
                            DE1(b,x) = DE1(b,x) + De1_type3{j,1,4,1,x}(param_s(:,m_now(j):m_next(j),:,:,t), param_s(:,m_now(j+1):m_next(j+1),:,:,t));
                        end
                    end
                    
                end
                
                if j == length(k1) - 1

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
                    param_s(1,b,1,1,t+1) = param_s(1,b,1,1,t) - eta_s1 * (DE1(1) + Ereg_coef * DEreg(1));
                    param_s(1,b,1,2,t+1) = param_s(1,b,1,2,t) - eta_s2 * (DE1(2) + Ereg_coef * DEreg(2));
                    param_s(1,b,1,3,t+1) = param_s(1,b,1,3,t) - eta_s3 * (DE1(3) + Ereg_coef * DEreg(3));

                end

            end
           
        end

        param_s(2,:,1,3,t+1) = param_s(1,:,1,3,t+1);
        param_s(1,:,2,3,t+1) = param_s(1,:,1,3,t+1);
        param_s(2,:,2,3,t+1) = param_s(1,:,1,3,t+1);

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