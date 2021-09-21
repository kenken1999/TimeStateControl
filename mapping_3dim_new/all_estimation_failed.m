clear;
close all;

tic

%---(l,m,n)=(0,0,0),(1,0,0),(0,1,0),(0,0,1)のfix, およびその他初期値の決定(線形補間)----------------------------------------------------

l_max = 3;
m_max = 3;
n_max = 3;

s = sym('s',[2 * l_max-1 2 * m_max-1 2 * n_max-1 3]); % l,m,nの順

iteration = 100;

param_s = zeros(2*l_max-1, 2*m_max-1, 2*n_max-1, 3, iteration);

param_s(1,1,1,:,1) = [1 1 pi/4];
param_s(2,1,1,:,1) = [1+1/sqrt(2) 1+1/sqrt(2) pi/4];
param_s(1,2,1,:,1) = [1 1 pi/2];
param_s(1,1,2,:,1) = [1-1/sqrt(2) 1+1/sqrt(2) pi/4];

s_l = param_s(2,1,1,:,1) - param_s(1,1,1,:,1);
s_m = param_s(1,2,1,:,1) - param_s(1,1,1,:,1);
s_n = param_s(1,1,2,:,1) - param_s(1,1,1,:,1);


for a = 1 : 2 * l_max - 1
    for b = 1 : 2 * m_max - 1
        for c = 1 : 2 * n_max - 1

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

dk1 = 1;   % 時間刻み
K1fin = 10;  %シミュレーション終了時間, length(k) = Kfin + 1
k1 = [0:dk1:K1fin];

u1_b1 = ones(length(k1),1) * 0.1; % 並進速度
u2_b1 = ones(length(k1),1) * (-pi/16); % 回転角速度

si_b1 = zeros(length(k1),3); % 観測するセンサ変数 , s = (s1, s2, s3) = (x ,y, θ)
si_b1(1,:) = [1 1 7*pi/12];    % (s1, s2, s3)の初期値を設定

si_c1 = zeros(length(k1),3); % 補正後のセンサ変数(zi,z3空間と等しい)、結果比較用
si_c1(1,:) = [0 0 pi/3];

l_now = zeros(length(k1),1);
m_now = zeros(length(k1),1);
n_now = zeros(length(k1),1);

l_next = zeros(length(k1),1);
m_next = zeros(length(k1),1);
n_next = zeros(length(k1),1);

l_real = zeros(length(k1),1);
m_real = zeros(length(k1),1);
n_real = zeros(length(k1),1);

rho_tmp = zeros(length(k1),3);
rho = zeros(length(k1),3);

break_switch = 0;


for a = 1 : 2*l_max - 1

    if break_switch == 1 
        break;
    end

    for b = 1 : 2*m_max - 1

        if break_switch == 1 
            break;
        end

        for c = 1 : 2*n_max - 1

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

            rho_tmp(1,:) = B \ x;

            if (0 <= rho_tmp(1,1)) && (rho_tmp(1,1) < 1)
                if (0 <= rho_tmp(1,2)) && (rho_tmp(1,2) < 1)
                    if (0 <= rho_tmp(1,3)) && (rho_tmp(1,3) < 1)

                        rho(1,:) = rho_tmp(1,:);

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


    for a = 1 : 2 * l_max - 1

        if break_switch == 1 
            break;
        end

        for b = 1 : 2 * m_max - 1

            if break_switch == 1 
                break;
            end

            for c = 1 : 2 * n_max - 1

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

                rho_tmp(j+1,:) = B \ x;

                if (0 <= rho_tmp(j+1,1)) && (rho_tmp(j+1,1) < 1)
                    if (0 <= rho_tmp(j+1,2)) && (rho_tmp(j+1,2) < 1)
                        if (0 <= rho_tmp(j+1,3)) && (rho_tmp(j+1,3) < 1)

                            rho(j+1,:) = rho_tmp(j+1,:);
    
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

    E1 = E1 + ( m_real(j) + P(2) - ((n_real(j+1) + P2(3)) - (n_real(j) + P(3))) / ((l_real(j+1) + P2(1)) - (l_real(j) + P(1))) ) ^ 2;

end



%---正則化項の追加-----------------------------

E4 = 0;

for a = 1 : l_max - 2
    for b = 1 : m_max - 2
        for c = 1 : n_max - 2

            E4 = E4 + (s(a+2,b,c,1) - 2 * s(a+1,b,c,1) - s(a,b,c,1)) ^ 2 + (s(a,b+2,c,1) - 2 * s(a,b+1,c,1) - s(a,b,c,1)) ^ 2 + (s(a,b,c+2,1) - 2 * s(a,b,c+1,1) - s(a,b,c,1)) ^ 2 ...
                    + (s(a+2,b,c,2) - 2 * s(a+1,b,c,2) - s(a,b,c,2)) ^ 2 + (s(a,b+2,c,2) - 2 * s(a,b+1,c,2) - s(a,b,c,2)) ^ 2 + (s(a,b,c+2,2) - 2 * s(a,b,c+1,2) - s(a,b,c,2)) ^ 2 ...
                    + (s(a+2,b,c,3) - 2 * s(a+1,b,c,3) - s(a,b,c,3)) ^ 2 + (s(a,b+2,c,3) - 2 * s(a,b+1,c,3) - s(a,b,c,3)) ^ 2 + (s(a,b,c+2,3) - 2 * s(a,b,c+1,3) - s(a,b,c,3)) ^ 2;

        end
    end
end

E1 = E1 + 0.00 * E4;



%---最急降下法によるパラメータの決定----------------------------

eta_s1 = 1.0 * 10 ^ (-8); % 学習率
eta_s2 = 1.0 * 10 ^ (-8);
eta_s3 = 1.0 * 10 ^ (-2);


%---Eの設定---------------------------

E1_initial = double(subs(E1, [s(:,:,:,:)],[param_s(:,:,:,:,1)]));
                             

disp('E1_initial = ')
disp(E1_initial)
disp('--------------------')


E1_value = zeros(1,iteration);


for a = 1 : 2 * l_max - 1
    for b = 1 : 2 * m_max - 1
        for c = 1 : 2 * n_max - 1

            DE1_s1(a,b,c) = diff(E1,s(a,b,c,1));
            DE1_s2(a,b,c) = diff(E1,s(a,b,c,2));
            DE1_s3(a,b,c) = diff(E1,s(a,b,c,3));

            DE1_s1_1{a,b,c} = matlabFunction(DE1_s1(a,b,c), 'vars', {s(:,:,:,:)});
            DE1_s2_1{a,b,c} = matlabFunction(DE1_s2(a,b,c), 'vars', {s(:,:,:,:)});
            DE1_s3_1{a,b,c} = matlabFunction(DE1_s3(a,b,c), 'vars', {s(:,:,:,:)});

            disp([a b c])
            disp('------------')

        end
    end
end



for t = 1:iteration - 1

    for a = 1 : 2 * l_max - 1
        for b = 1 : 2 * m_max - 1
            for c = 1 : 2 * n_max - 1

            DE1_s1_2 = DE1_s1_1{a,b,c}(param_s(:,:,:,:,t));
            DE1_s2_2 = DE1_s2_1{a,b,c}(param_s(:,:,:,:,t));
            DE1_s3_2 = DE1_s3_1{a,b,c}(param_s(:,:,:,:,t));

            param_s(a,b,c,1,t+1) = param_s(a,b,c,1,t) - eta_s1 * DE1_s1_2;
            param_s(a,b,c,2,t+1) = param_s(a,b,c,2,t) - eta_s2 * DE1_s2_2;
            param_s(a,b,c,3,t+1) = param_s(a,b,c,3,t) - eta_s3 * DE1_s3_2;

            end
        end
    end


    E1_value(t) = double(subs(E1, [s(:,:,:,:)],[param_s(:,:,:,:,t+1)]));


    disp('t = ')
    disp(t)
    disp('E1(t) = ')
    disp(E1_value(t))
 

    if t > 1

        if (E1_value(t) > E1_value(t-1))
            disp('Ef1_2が増加しました')
            disp('iterationを強制終了します')
            iteration = t-1;
            break;
        end
        if t == iteration - 1
            iteration = t+1;
            disp('iterationを正常に終了することができました！')
            break
        end

    end

    disp('--------------------')

end

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

rho_2 = zeros(length(k1),3);


for j = 1:length(k1)

    for a = 1 : 2 * l_max - 1

        if break_switch == 1 
            break;
        end

        for b = 1 : 2 * m_max - 1

            if break_switch == 1 
                break;
            end

            for c = 1 : 2 * n_max - 1

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

                rho_tmp(j,:) = B \ x;

                if (0 <= rho_tmp(j,1)) && (rho_tmp(j,1) < 1)
                    if (0 <= rho_tmp(j,2)) && (rho_tmp(j,2) < 1)
                        if (0 <= rho_tmp(j,3)) && (rho_tmp(j,3) < 1)

                            rho_2(j,:) = rho_tmp(j,:);
    
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

    break_switch = 0;

    z1_b1(j) = l_real_2(j) + rho_2(j,1);
    z2_b1(j) = m_real_2(j) + rho_2(j,2);
    z3_b1(j) = n_real_2(j) + rho_2(j,3);

end


% 推定結果のplot--------------------------------------

figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c1(:,3), si_c1(:,1), '--m', si_c1(:,3), z1_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z1 = f1(s) = s1 の答え合わせ
xlabel("s3' = θ")
ylabel('z1')
legend("真値：s1'",'推定値：z1')

hold off;

figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c1(:,3), tan(si_c1(:,3)), '--m', si_c1(:,3), z2_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z1 = f1(s) = s1 の答え合わせ
xlabel("s3' = θ")
ylabel("z2")
legend("真値：tan(s3')",'推定値：z2')

hold off;

figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c1(:,3), si_c1(:,2), '--m', si_c1(:,3), z3_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z3 = f3(s) = s2 の答え合わせ
xlabel("s3' = θ")
ylabel('z3')
legend("真値：s2'",'推定値：z3')

hold off;



%---p,qの初期値の決定--------------------------------------------------------

p_max = 3;
q_max = 3;

sp = sym('sp',[p_max 3]);
sq = sym('sq',[q_max 3]);

iteration2 = 100;

param_sp = zeros(p_max, 3, iteration2);
param_sq = zeros(q_max, 3, iteration2);


for p = 1 : p_max
    param_sp(p,:,1) = param_sp(1,:,1) + (p-1) * [0.2 0.2 0.2];
end

for q = 1 : q_max
    param_sq(q,:,1) = param_sq(1,:,1) + (q-1) * [0.2 0.2 0.2];
end


%---誤差関数の定義(g,h)----------------------------------------------------------

p_now = zeros(length(k1),1);
q_now = zeros(length(k1),1);

p_next = zeros(length(k1),1);
q_next = zeros(length(k1),1);

p_real = zeros(length(k1),1);
q_real = zeros(length(k1),1);

rho_g_tmp = zeros(length(k1),1);
rho_g = zeros(length(k1),1);

rho_h_tmp = zeros(length(k1),1);
rho_h = zeros(length(k1),1);

break_switch = 0;


for p = 1 : p_max - 1

    if break_switch == 1 
        break;
    end

    p2 = p + 1;
    p_real(1) = p - 1;

    rho_g_tmp(1) = ((si_b1(1,1) - param_sp(p,1,1)) / (param_sp(p2,1,1) - param_sp(p,1,1)))...
                 + ((si_b1(1,2) - param_sp(p,2,1)) / (param_sp(p2,2,1) - param_sp(p,2,1)))...
                 + ((si_b1(1,3) - param_sp(p,3,1)) / (param_sp(p2,3,1) - param_sp(p,3,1)));

    if (0 <= rho_g_tmp(1)) && (rho_g_tmp(1) > 1)

        rho_g(1) = rho_g_tmp(1);

        p_now(1) = p;
        p_next(1) = p2;

        break_switch = 1;

    end

end

break_switch = 0;

for q = 1 : q_max - 1

    if break_switch == 1 
        break;
    end

    q2 = q + 1;
    q_real(1) = q - 1;

    rho_h_tmp(1) = ((si_b1(1,1) - param_sq(q,1,1)) / (param_sq(q2,1,1) - param_sq(q,1,1)))...
                 + ((si_b1(1,2) - param_sq(q,2,1)) / (param_sq(q2,2,1) - param_sq(q,2,1)))...
                 + ((si_b1(1,3) - param_sq(q,3,1)) / (param_sq(q2,3,1) - param_sq(q,3,1)));

    if (0 <= rho_h_tmp(1)) && (rho_h_tmp(1) > 1)

        rho_h(1) = rho_h_tmp(1);

        q_now(1) = q;
        q_next(1) = q2;

        break_switch = 1;

    end

end

break_switch = 0;

E2 = 0;
E3 = 0;


for j = 1 : length(k1) - 1
   

    for p = 1 : p_max - 1

        if break_switch == 1 
            break;
        end
    
        p2 = p + 1;
        p_real(j+1) = p - 1;
    
        rho_g_tmp(j+1) = ((si_b1(j+1,1) - param_sp(p,1,1)) / (param_sp(p2,1,1) - param_sp(p,1,1)))...
                       + ((si_b1(j+1,2) - param_sp(p,2,1)) / (param_sp(p2,2,1) - param_sp(p,2,1)))...
                       + ((si_b1(j+1,3) - param_sp(p,3,1)) / (param_sp(p2,3,1) - param_sp(p,3,1)));
    
        if (0 <= rho_g_tmp(j+1)) && (rho_g_tmp(j+1) > 1)
    
            rho_g(j+1) = rho_g_tmp(j+1);
    
            p_now(j+1) = p;
            p_next(j+1) = p2;
    
            break_switch = 1;
    
        end
    
    end

    break_switch = 0;

    for q = 1 : q_max - 1

        if break_switch == 1 
            break;
        end
    
        q2 = q + 1;
        q_real(j+1) = q - 1;
    
        rho_h_tmp(j+1) = ((si_b1(j+1,1) - param_sq(q,1,1)) / (param_sq(q2,1,1) - param_sq(q,1,1)))...
                       + ((si_b1(j+1,2) - param_sq(q,2,1)) / (param_sq(q2,2,1) - param_sq(q,2,1)))...
                       + ((si_b1(j+1,3) - param_sq(q,3,1)) / (param_sq(q2,3,1) - param_sq(q,3,1)));
    
        if (0 <= rho_h_tmp(j+1)) && (rho_h_tmp(j+1) > 1)
    
            rho_h(j+1) = rho_h_tmp(j+1);
    
            q_now(j+1) = q;
            q_next(j+1) = q2;
    
            break_switch = 1;
    
        end
    
    end
    
    break_switch = 0;

    p = p_now(j);
    q = q_now(j);

    p_p = p_next(j);
    q_p = q_next(j);

    p2 = p_now(j+1);
    q2 = q_now(j+1);
 
    p2_p = p_next(j+1);
    q2_p = q_next(j+1);


    Pg = ((si_b1(j,1) - sp(p,1)) / (sp(p_p,1) - sp(p,1)))...
       + ((si_b1(j,2) - sp(p,2)) / (sp(p_p,2) - sp(p,2)))...
       + ((si_b1(j,3) - sp(p,3)) / (sp(p_p,3) - sp(p,3)));

    Pg2 = ((si_b1(j+1,1) - sp(p2,1)) / (sp(p2_p,1) - sp(p2,1)))...
        + ((si_b1(j+1,2) - sp(p2,2)) / (sp(p2_p,2) - sp(p2,2)))...
        + ((si_b1(j+1,3) - sp(p2,3)) / (sp(p2_p,3) - sp(p2,3)));


    Ph = ((si_b1(j,1) - sq(q,1)) / (sq(q_p,1) - sq(q,1)))...
       + ((si_b1(j,2) - sq(q,2)) / (sq(q_p,2) - sq(q,2)))...
       + ((si_b1(j,3) - sq(q,3)) / (sq(q_p,3) - sq(q,3)));

    Ph2 = ((si_b1(j+1,1) - sq(q2,1)) / (sq(q2_p,1) - sq(q2,1)))...
        + ((si_b1(j+1,2) - sq(q2,2)) / (sq(q2_p,2) - sq(q2,2)))...
        + ((si_b1(j+1,3) - sq(q2,3)) / (sq(q2_p,3) - sq(q2,3)));


    E2 = E2 + ( u2_b1(j) / u1_b1(j) * (p_real(j) + Pg) - ((z2_b1(j+1) - z2_b1(j)) / (z1_b1(j+1) - z1_b1(j))) ) ^ 2;

    E3 = E3 + ( u1_b1(j) * (q_real(j) + Ph) - ((z1_b1(j+1) - z1_b1(j)) / dk1) ) ^ 2;

end


%---最急降下法によるパラメータの決定----------------------------

eta_sp1 = 1.0 * 10 ^ (-1); % 学習率
eta_sp2 = 1.0 * 10 ^ (-1);
eta_sp3 = 1.0 * 10 ^ (-1);

eta_sq1 = 1.0 * 10 ^ (-1); % 学習率
eta_sq2 = 1.0 * 10 ^ (-1);
eta_sq3 = 1.0 * 10 ^ (-1);


%---Eの設定---------------------------

E2_initial = double(subs(E2, [sp(:,:)],[param_sp(:,:,1)]));
E3_initial = double(subs(E3, [sq(:,:)],[param_sq(:,:,1)]));
                             

disp('E2_initial = ')
disp(E2_initial)
disp('E3_initial = ')
disp(E3_initial)
disp('--------------------')


E2_value = zeros(1,iteration);
E3_value = zeros(1,iteration);


for p = 1 : p_max

    DE2_sp1(p) = diff(E2,sp(p,1));
    DE2_sp2(p) = diff(E2,sp(p,2));
    DE2_sp3(p) = diff(E2,sp(p,3));

    DE2_sp1_1{p} = matlabFunction(DE2_sp1(p), 'vars', {sp(:,:)});
    DE2_sp2_1{p} = matlabFunction(DE2_sp2(p), 'vars', {sp(:,:)});
    DE2_sp3_1{p} = matlabFunction(DE2_sp3(p), 'vars', {sp(:,:)});

    disp(p)
    disp('------------')

end

for q = 1 : q_max

    DE3_sq1(q) = diff(E3,sq(q,1));
    DE3_sq2(q) = diff(E3,sq(q,2));
    DE3_sq3(q) = diff(E3,sq(q,3));

    DE3_sq1_1{q} = matlabFunction(DE3_sq1(q), 'vars', {sq(:,:)});
    DE3_sq2_1{q} = matlabFunction(DE3_sq2(q), 'vars', {sq(:,:)});
    DE3_sq3_1{q} = matlabFunction(DE3_sq3(q), 'vars', {sq(:,:)});

    disp(q)
    disp('------------')

end



for t = 1:iteration2 - 1

    for p = 1 : p_max

        DE2_sp1_2 = DE2_sp1_1{p}(param_sp(:,:,t));
        DE2_sp2_2 = DE2_sp2_1{p}(param_sp(:,:,t));
        DE2_sp3_2 = DE2_sp3_1{p}(param_sp(:,:,t));

        param_sp(p,1,t+1) = param_sp(p,1,t) - eta_sp1 * DE2_sp1_2;
        param_sp(p,2,t+1) = param_sp(p,2,t) - eta_sp2 * DE2_sp2_2;
        param_sp(p,3,t+1) = param_sp(p,3,t) - eta_sp3 * DE2_sp3_2;
          
    end

    for q = 1 : q_max

        DE3_sq1_2 = DE3_sq1_1{q}(param_sq(:,:,t));
        DE3_sq2_2 = DE3_sq2_1{q}(param_sq(:,:,t));
        DE3_sq3_2 = DE3_sq3_1{q}(param_sq(:,:,t));

        param_sq(q,1,t+1) = param_sq(q,1,t) - eta_sq1 * DE3_sq1_2;
        param_sq(q,2,t+1) = param_sq(q,2,t) - eta_sq2 * DE3_sq2_2;
        param_sq(q,3,t+1) = param_sq(q,3,t) - eta_sq3 * DE3_sq3_2;
          
    end


    E2_value(t) = double(subs(E2, [sp(:,:)],[param_sp(:,:,t+1)]));
    E3_value(t) = double(subs(E3, [sq(:,:)],[param_sq(:,:,t+1)]));


    disp('t = ')
    disp(t)
    disp('E2(t) = ')
    disp(E2_value(t))
    disp('E3(t) = ')
    disp(E3_value(t))
 

    if t > 1

        if (E2_value(t) > E2_value(t-1))
            disp('E2が増加しました')
            disp('iteration2を強制終了します')
            iteration2 = t-1;
            break;
        end
        if (E3_value(t) > E3_value(t-1))
            disp('E3が増加しました')
            disp('iteration2を強制終了します')
            iteration2 = t-1;
            break;
        end
        if t == iteration2 - 1
            iteration2 = t+1;
            disp('iteration2を正常に終了することができました！')
            break
        end

    end

    disp('--------------------')

end


% z1,z2,z3の推定結果取得--------------------------------------

g_b1 = zeros(length(k1),1);
h_b1 = zeros(length(k1),1);

p_now_2 = zeros(length(k1),1);
q_now_2 = zeros(length(k1),1);

p_next_2 = zeros(length(k1),1);
q_next_2 = zeros(length(k1),1);

p_real_2 = zeros(length(k1),1);
q_real_2 = zeros(length(k1),1);

rho_g_2 = zeros(length(k1),1);
rho_h_2 = zeros(length(k1),1);


for j = 1:length(k1)

    for p = 1 : p_max - 1

        if break_switch == 1 
            break;
        end

        p2 = p + 1;
        p_real(j+1) = p - 1;

        rho_g_tmp(j) = ((si_b1(j,1) - param_sp(p,1,iteration2)) / (param_sp(p2,1,iteration2) - param_sp(p,1,iteration2)))...
                     + ((si_b1(j,2) - param_sp(p,2,iteration2)) / (param_sp(p2,2,iteration2) - param_sp(p,2,iteration2)))...
                     + ((si_b1(j,3) - param_sp(p,3,iteration2)) / (param_sp(p2,3,iteration2) - param_sp(p,3,iteration2)));

        if (0 <= rho_g_tmp(j)) && (rho_g_tmp(j) < 1)

            rho_g_2(j) = rho_g_tmp(j);

            p_now_2(j) = p;
            p_next_2(j) = p2;

            break_switch = 1;

        end

    end

    break_switch = 0;


    for q = 1 : q_max - 1

        if break_switch == 1 
            break;
        end

        q2 = q + 1;
        q_real(j+1) = q - 1;

        rho_h_tmp(j) = ((si_b1(j,1) - param_sq(q,1,iteration2)) / (param_sq(q2,1,iteration2) - param_sq(q,1,iteration2)))...
                     + ((si_b1(j,2) - param_sq(q,2,iteration2)) / (param_sq(q2,2,iteration2) - param_sq(q,2,iteration2)))...
                     + ((si_b1(j,3) - param_sq(q,3,iteration2)) / (param_sq(q2,3,iteration2) - param_sq(q,3,iteration2)));

        if (0 <= rho_h_tmp(j)) && (rho_h_tmp(j) < 1)

            rho_h_2(j) = rho_h_tmp(j);

            q_now_2(j) = q;
            q_next_2(j) = q2;

            break_switch = 1;

        end

    end

    break_switch = 0;

    g_b1(j) = p_real_2(j) + rho_g_2(j);
    h_b1(j) = q_real_2(j) + rho_h_2(j);

end


toc

% 推定結果のplot--------------------------------------

figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

for i = 1 : length(k1)
    g_ans(i) = 1 / (cos(si_c1(i,3)) * cos(si_c1(i,3)) * cos(si_c1(i,3)));
end

plot(si_c1(:,3), g_ans(:), '--m', si_c1(:,3), g_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("s3' = θ")
ylabel('g')
legend("真値：1/cos^3(s3')",'推定値：g')

hold off;


figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c1(:,3), cos(si_c1(:,3)), '--m', si_c1(:,3), h_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel("s3' = θ")
ylabel("h")
legend("真値：cos(s3')",'推定値：h')

hold off;
