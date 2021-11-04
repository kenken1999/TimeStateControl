clear;
close all;

tic

%---(l,m,n)=(0,0,0),(1,0,0),(0,1,0),(0,0,1)のfix, およびその他初期値の決定(線形補間)----------------------------------------------------

l_max = 2;
m_max = 10;
n_max = 2;

iteration = 100;

s = sym('s',[l_max m_max n_max 3]); % l,m,nの順

param_s = zeros(l_max, m_max, n_max, 3, iteration);

param_s(1,1,1,:,1) = [1 1 pi/4];
param_s(2,1,1,:,1) = [1+1/sqrt(2) 1+1/sqrt(2) pi/4];
param_s(1,2,1,:,1) = [1 1 pi/3];
param_s(1,1,2,:,1) = [1-1/sqrt(2) 1+1/sqrt(2) pi/4];

s_l = param_s(2,1,1,:,1) - param_s(1,1,1,:,1);
s_m = param_s(1,2,1,:,1) - param_s(1,1,1,:,1);
s_n = param_s(1,1,2,:,1) - param_s(1,1,1,:,1);

l_max_now = 2;
m_max_now = 3;
n_max_now = 2;

m_start_change = 1;
m_save = 1;

max_switch = 0;


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


%---サンプル収集と誤差関数の定義----------------------------------------------------

imax = 14;

for i = 1 : imax

    dk1 = 0.1;   % 時間刻み
    K1fin = 1.9;  %シミュレーション終了時間, length(k) = Kfin + 1
    k1 = [0:dk1:K1fin];

    u1_b1 = ones(length(k1),1) * 0.5; % 並進速度

    if rem(i,2) == 1
        u2_b1 = ones(length(k1),1) * (0.6); % 回転角速度
    else
        u2_b1 = ones(length(k1),1) * (0.6); % 回転角速度
    end

    si_b1 = zeros(length(k1),3); % 観測するセンサ変数 , s = (s1, s2, s3) = (x ,y, θ)
    si_b1(1,:) = [1 1 pi/4];    % (s1, s2, s3)の初期値を設定

    si_c1 = zeros(length(k1),3); % 補正後のセンサ変数(zi,z3空間と等しい)、結果比較用
    si_c1(1,:) = [0 0 0];

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

        if rem(i,2) == 1 && m_max_now < m_max

            s_m = param_s(1,m_max_now,1,:,1) - param_s(1,m_max_now - 1,1,:,1);

            for b = m_max_now : m_max

                if b <= m_max
                    m_coef = b - m_max_now;
                end
    
                param_s(1,b,1,:,1) = param_s(1,m_max_now,1,:,1) + m_coef * s_m;
                              
            end

        end

        param_s(2,:,1,3,1) = param_s(1,:,1,3,1);
        param_s(1,:,2,3,1) = param_s(1,:,1,3,1);
        param_s(2,:,2,3,1) = param_s(1,:,1,3,1);

    end

    break_switch = 0;


    for a = 1 : l_max

        if break_switch == 1 
            break;
        end

        for b = 1 : m_max

            if break_switch == 1 
                break;
            end

            for c = 1 : n_max

                if break_switch == 1 
                    break;
                end

                if (a == l_max) || (b == m_max) || (c == n_max)
                    continue;
                end

                if a < l_max
                    a2 = a + 1;
                    l_real(1) = a - 1;
                elseif a == l_max + 1
                    a2 = 1;
                    l_real(1) = -1;
                else
                    a2 = a - 1;
                    l_real(1) = - a + l_max;
                end

                if b < m_max
                    b2 = b + 1;
                    m_real(1) = b - 1;
                elseif b == m_max + 1
                    b2 = 1;
                    m_real(1) = -1;
                else
                    b2 = b - 1;
                    m_real(1) = - b + m_max;
                end

                if c < n_max
                    c2 = c + 1;
                    n_real(1) = c - 1;
                elseif c == n_max + 1
                    c2 = 1;
                    n_real(1) = -1;
                else
                    c2 = c - 1;
                    n_real(1) = - c + n_max;
                end
            

                A = [param_s(a2,b,c,:,1) - param_s(a,b,c,:,1); param_s(a,b2,c,:,1) - param_s(a,b,c,:,1); param_s(a,b,c2,:,1) - param_s(a,b,c,:,1)];

                B = transpose(reshape(A,[3,3]));

                x = [si_b1(1,1) - param_s(a,b,c,1,1); si_b1(1,2) - param_s(a,b,c,2,1); si_b1(1,3) - param_s(a,b,c,3,1)];

                rho_tmp(1,a,b,c,:) = B \ x;

                if (0 <= rho_tmp(1,a,b,c,1)) && (rho_tmp(1,a,b,c,1) <= 1)
                    if (0 <= rho_tmp(1,a,b,c,2)) && (rho_tmp(1,a,b,c,2) <= 1)
                        if (0 <= rho_tmp(1,a,b,c,3)) && (rho_tmp(1,a,b,c,3) <= 1)

                            rho(1,:) = rho_tmp(1,a,b,c,:);

                            l_now(1) = a;
                            m_now(1) = b;
                            n_now(1) = c;

                            l_next(1) = a2;
                            m_next(1) = b2;
                            n_next(1) = c2;

                            break_switch = 1;

                        end
                    end
                end
                
            end

        end

    end

    break_switch = 0;

    E1 = 0;


    for j = 1 : length(k1) - 1
    
        si_b1(j+1,3) = si_b1(j,3) + u2_b1(j+1) * dk1;
        si_b1(j+1,1) = si_b1(j,1) + u1_b1(j+1) * cos(si_b1(j+1,3)) * dk1;
        si_b1(j+1,2) = si_b1(j,2) + u1_b1(j+1) * sin(si_b1(j+1,3)) * dk1;

        si_c1(j+1,3) = si_c1(j,3) + u2_b1(j+1) * dk1;
        si_c1(j+1,1) = si_c1(j,1) + u1_b1(j+1) * cos(si_c1(j+1,3)) * dk1;
        si_c1(j+1,2) = si_c1(j,2) + u1_b1(j+1) * sin(si_c1(j+1,3)) * dk1;


        for a = 1 : l_max

            if break_switch == 1 
                break;
            end

            for b = 1 : m_max

                if break_switch == 1 
                    break;
                end

                for c = 1 : n_max

                    if break_switch == 1 
                        break;
                    end

                    if (a == l_max) || (b == m_max) || (c == n_max)
                        continue;
                    end
        
                    if a < l_max
                        a2 = a + 1;
                        l_real(j+1) = a - 1;
                    elseif a == l_max + 1
                        a2 = 1;
                        l_real(j+1) = -1;
                    else
                        a2 = a - 1;
                        l_real(j+1) = - a + l_max; 
                    end
        
                    if b < m_max
                        b2 = b + 1;
                        m_real(j+1) = b - 1;
                    elseif b == m_max + 1
                        b2 = 1;
                        m_real(j+1) = -1;
                    else
                        b2 = b - 1;
                        m_real(j+1) = - b + m_max;
                    end
        
                    if c < n_max
                        c2 = c + 1;
                        n_real(j+1) = c - 1;
                    elseif c == n_max + 1
                        c2 = 1;
                        n_real(j+1) = -1;
                    else
                        c2 = c - 1;
                        n_real(j+1) = - c + n_max;
                    end
        
                    A = [param_s(a2,b,c,:,1) - param_s(a,b,c,:,1); param_s(a,b2,c,:,1) - param_s(a,b,c,:,1); param_s(a,b,c2,:,1) - param_s(a,b,c,:,1)];

                    B = transpose(reshape(A,[3,3]));

                    x = [si_b1(j+1,1) - param_s(a,b,c,1,1); si_b1(j+1,2) - param_s(a,b,c,2,1); si_b1(j+1,3) - param_s(a,b,c,3,1)];

                    rho_tmp(j+1,a,b,c,:) = B \ x;

                    if (0 <= rho_tmp(j+1,a,b,c,1)) && (rho_tmp(j+1,a,b,c,1) <= 1)
                        if (0 <= rho_tmp(j+1,a,b,c,2)) && (rho_tmp(j+1,a,b,c,2) <= 1)
                            if (0 <= rho_tmp(j+1,a,b,c,3)) && (rho_tmp(j+1,a,b,c,3) <= 1)

                                rho(j+1,:) = rho_tmp(j+1,a,b,c,:);
        
                                l_now(j+1) = a;
                                m_now(j+1) = b;
                                n_now(j+1) = c;

                                l_next(j+1) = a2;
                                m_next(j+1) = b2;
                                n_next(j+1) = c2;

                                break_switch = 1;
        
                            end
                        end
                    end

                end
            end
        end

        %欠損時の補間
        if (break_switch == 0)

            l_now(j+1) = l_now(j);
            m_now(j+1) = m_now(j);
            n_now(j+1) = n_now(j);

            l_next(j+1) = l_next(j);
            m_next(j+1) = m_next(j);
            n_next(j+1) = n_next(j);

            l_real(j+1) = l_real(j);
            m_real(j+1) = m_real(j);
            n_real(j+1) = n_real(j);

            a = l_now(j+1);
            b = m_now(j+1);
            c = n_now(j+1);

            rho(j+1,:) = rho_tmp(j+1,a,b,c,:);
            
            disp(j)
            disp("補間しました")

        end

        break_switch = 0;

        l = l_now(j);
        m = m_now(j);
        n = n_now(j);

        l_p = l_next(j);
        m_p = m_next(j);
        n_p = n_next(j);

        l2 = l_now(j+1);
        m2 = m_now(j+1);
        n2 = n_now(j+1);

        l2_p = l_next(j+1);
        m2_p = m_next(j+1);
        n2_p = n_next(j+1);


        G = [s(l_p,m,n,:) - s(l,m,n,:); s(l,m_p,n,:) - s(l,m,n,:); s(l,m,n_p,:) - s(l,m,n,:)];
        H =  transpose(reshape(G,[3,3]));
        y = [si_b1(j,1) - s(l,m,n,1); si_b1(j,2) - s(l,m,n,2); si_b1(j,3) - s(l,m,n,3)];
        P = H \ y;

        G2 = [s(l2_p,m2,n2,:) - s(l2,m2,n2,:); s(l2,m2_p,n2,:) - s(l2,m2,n2,:); s(l2,m2,n2_p,:) - s(l2,m2,n2,:)];
        H2 =  transpose(reshape(G2,[3,3]));
        y2 = [si_b1(j+1,1) - s(l2,m2,n2,1); si_b1(j+1,2) - s(l2,m2,n2,2); si_b1(j+1,3) - s(l2,m2,n2,3)];
        P2 = H2 \ y2;

        E1 = E1 + ( tan(pi/12) * (m_real(j) + P(2)) - ((n_real(j+1) + P2(3)) - (n_real(j) + P(3))) / ((l_real(j+1) + P2(1)) - (l_real(j) + P(1))) ) ^ 2;


    end

    disp('-----')



    %---正則化項の追加-----------------------------

    if rem(i,2) == 1
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



    E4 = 0;

    for a = 1 : 2
        for b = 1 : m_max - 2
            for c = 1 : 2

                E4 = E4 + (( (s(a,b+2,c,1) - s(a,b+1,c,1)) ^ 2 + (s(a,b+2,c,2) - s(a,b+1,c,2)) ^ 2 + (s(a,b+2,c,3) - s(a,b+1,c,3)) ^ 2 )... 
                        - ( (s(a,b+1,c,1) - s(a,b,c,1)) ^ 2 + (s(a,b+1,c,2) - s(a,b,c,2)) ^ 2 + (s(a,b+1,c,3) - s(a,b,c,3)) ^ 2 )) ^ 2;
                    
            end
        end
    end


    E1_initial = double(subs(E1, [s(:,:,:,:)],[param_s(:,:,:,:,1)]));    
    E4_initial = double(subs(E4, [s(:,:,:,:)],[param_s(:,:,:,:,1)]));

    E4_coef = 10;

    disp('E4_initial = ')
    disp(E4_initial)
    disp('E4_coef = ')
    disp(E4_coef)

    E1 = E1 + E4_coef * E4;


    %---最急降下法のパラメータ決定----------------------------

    eta_s1 = 0.0 * 10 ^ (-6); % 学習率
    eta_s2 = 0.0 * 10 ^ (-6);
    eta_s3 = 1.0 * 10 ^ (-2);

    iteration = 100;

    stop_switch = 0;

    %---Eの設定---------------------------

    E1_initial = double(subs(E1, [s(:,:,:,:)],[param_s(:,:,:,:,1)]));                 

    disp('E1_initial = ')
    disp(E1_initial)
    disp('--------------------')


    E1_value = zeros(1,iteration);

    m_start_change = m_save;

    if m_max_now <= m_max && max_switch < 2
        m_start_change = 1;
        m_save = m_start_change;
        m_max_change = m_max_now;
        if m_max_now == m_max
            max_switch = max_switch + 1;
        end
    else
        if rem(i,2) == 1
            m_start_change = m_start_change + 1;
            m_save = m_start_change;
            m_max_change = m_max;
        end
    end


    DE1_s1_1 = cell(l_max, m_max, n_max);
    DE1_s2_1 = cell(l_max, m_max, n_max);
    DE1_s3_1 = cell(l_max, m_max, n_max);


    for a = l_max_now - 1 : l_max_now
        for b = m_start_change : m_max_change
            for c = n_max_now - 1 : n_max_now

                % DE1_s1(a,b,c) = diff(E1,s(a,b,c,1));
                % DE1_s2(a,b,c) = diff(E1,s(a,b,c,2));
                % DE1_s3(a,b,c) = diff(E1,s(a,b,c,3));

                % DE1_s1_1{a,b,c} = matlabFunction(DE1_s1(a,b,c), 'vars', {s(:,:,:,:)});
                % DE1_s2_1{a,b,c} = matlabFunction(DE1_s2(a,b,c), 'vars', {s(:,:,:,:)});
                % DE1_s3_1{a,b,c} = matlabFunction(DE1_s3(a,b,c), 'vars', {s(:,:,:,:)});

                DE1_s1_1{a,b,c} = matlabFunction(diff(E1,s(a,b,c,1)), 'vars', {s(:,:,:,:)});
                DE1_s2_1{a,b,c} = matlabFunction(diff(E1,s(a,b,c,2)), 'vars', {s(:,:,:,:)});
                DE1_s3_1{a,b,c} = matlabFunction(diff(E1,s(a,b,c,3)), 'vars', {s(:,:,:,:)});

                if a == 1 && b == m_start_change && c == 1
                    disp('i = ')
                    disp(i)
                    disp('m_start_change = ')
                    disp(m_start_change)
                    disp('m_max_change = ')
                    disp(m_max_change)
                    disp([a b c])
                    disp('-----')
                end

            end
        end
    end



    for t = 1:iteration - 1

        param_s(:,:,:,:,t+1) = param_s(:,:,:,:,t);

        for a = l_max_now - 1 : l_max_now
            for b = m_start_change : m_max_change
                for c = n_max_now - 1 : n_max_now

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


        E1_value(t) = double(subs(E1, (s(:,:,:,:)),(param_s(:,:,:,:,t+1))));

        if (t == 1 || t == 50 || t == 99)

            disp('i = ')
            disp(i)
            disp('t = ')
            disp(t)
            disp('E1(t) = ')
            disp(E1_value(t))

            disp('--------------------')

        end

        if t > 1

            if (E1_value(t) > E1_value(t-1))
                disp('t = ')
                disp(t)
                disp('Ef1が増加しました')
                disp('学習率を下げて再開します')
                disp('--------------------')
                eta_s1 = eta_s1 * 0.5; 
                eta_s2 = eta_s2 * 0.5;
                eta_s3 = eta_s3 * 0.5;
                param_s(:,:,:,:,t+1) = param_s(:,:,:,:,t);
                stop_switch = stop_switch + 1;
            end
            if stop_switch == 8
                disp('t = ')
                disp(t)
                disp('Ef1が増加しました')
                disp('iterationを強制終了します')
                disp("####################")
                iteration = t-1;
                break;
            end
            if t == iteration - 1
                iteration = t+1;
                disp('iterationを正常に終了することができました！')
                disp("####################")
                break
            end

        end
    
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

            for a = 1 : l_max

                if break_switch == 1 
                    break;
                end

                for b = 1 : m_max

                    if break_switch == 1 
                        break;
                    end

                    for c = 1 : n_max

                        if break_switch == 1 
                            break;
                        end

                        if (a == l_max) || (b == m_max) || (c == n_max)
                            continue;
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