clear;
close all;

tic

%z1&z3軸の推定(k1軸上)----------------------------------------------------

dk1 = 1;   %時間刻み
K1fin = 14; %シミュレーション終了時間, length(k) = Kfin + 1
k1 = [0:dk1:K1fin];

u1_b1 = zeros(length(k1),1); %並進速度
u2_b1 = zeros(length(k1),1); %回転角速度

si_b1 = zeros(length(k1),3); %観測するセンサ変数 , s = (s1, s2, s3) = (x ,y, θ)
si_b1(1,:) = [1 1 -pi/4];    %(s1, s2, s3)の初期値を設定

si_c1 = zeros(length(k1),3); %補正後のセンサ変数(zi,z3空間と等しい)、結果比較用
si_c1(1,:) = [0 -1 -pi/2];


f1_b1 = zeros(length(k1),1);   %写像f1:s→z1の推定
f3_b1 = zeros(length(k1),1);   %写像f3:s→z3の推定


l1_now = zeros(length(k1),1);  % l_now = 時刻kのl, 値が飛び飛び or 被る可能性あり, 0や負の値になる可能性あり
m1_now = zeros(length(k1),1);  % m_now = 時刻kのl, 値が飛び飛び or 被る可能性あり, 0や負の値になる可能性あり
n1_now = zeros(length(k1),1);  % n_now = 時刻kのl, 値が飛び飛び or 被る可能性あり, 0や負の値になる可能性あり

l1_num = ones(length(k1),1); %推定した格子点に前から順に1,2,3,...と番号をつけていく
m1_num = ones(length(k1),1);
n1_num = ones(length(k1),1);

rho1_1 = zeros(length(k1),1);
rho2_1 = zeros(length(k1),1);
rho3_1 = zeros(length(k1),1);

alpha = sym('alpha',[10 10 10]); %l,m,nの順
gamma = sym('gamma',[10 10 10]);

sigma1 = 0.2; %s1のスケーリング定数
sigma2 = 0.2; %s1のスケーリング定数
sigma3 = pi/4; %s1のスケーリング定数

Ef1 = 0;
Ef3 = 0;
E4_f1 = 0;
E4_f3 = 0;

l_max = 1; %現在の番号付けの最大値を保存， 全軸共通
m_max = 1;
n_max = 1;


l1_now(1) = floor(si_b1(1,1) / sigma1);
rho1_1(1) = si_b1(1,1) / sigma1  - l1_now(1);

m1_now(1) = floor(si_b1(1,2) / sigma2);
rho2_1(1) = si_b1(1,2) / sigma2  - m1_now(1);

n1_now(1) = floor(si_b1(1,3) / sigma3);
rho3_1(1) = si_b1(1,3) / sigma3  - n1_now(1);


y1 = zeros(length(k1),1);
y1(1) = -1;


%---誤差関数の定義（k1軸上）----------------------------------------------------


for j = 1 : length(k1) - 1

    if mod(j,10) == 0
        u1_b1(j+1) = -1;
        u2_b1(j+1) = 0;
        y1(j+1) = y1(j) + 1;
    elseif mod(j,5) == 0
        u1_b1(j+1) = 1;
        u2_b1(j+1) = 0;
        y1(j+1) = y1(j) + 1;
    elseif j < 5
        u1_b1(j+1) = 0;
        u2_b1(j+1) = pi/4;
        y1(j+1) = y1(j);
    elseif j > 5 && j < 10 
        u1_b1(j+1) = 0;
        u2_b1(j+1) = -pi/4;
        y1(j+1) = y1(j);
    else
        u1_b1(j+1) = 0;
        u2_b1(j+1) = pi/4;
        y1(j+1) = y1(j);
    end

   
    si_b1(j+1,3) = si_b1(j,3) + u2_b1(j+1) * dk1;
    si_b1(j+1,1) = si_b1(j,1) + u1_b1(j+1) * cos(si_b1(j+1,3)) * dk1;
    si_b1(j+1,2) = si_b1(j,2) + u1_b1(j+1) * sin(si_b1(j+1,3)) * dk1;

    si_c1(j+1,3) = si_c1(j,3) + u2_b1(j+1) * dk1;
    si_c1(j+1,1) = si_c1(j,1) + u1_b1(j+1) * cos(si_c1(j+1,3)) * dk1;
    si_c1(j+1,2) = si_c1(j,2) + u1_b1(j+1) * sin(si_c1(j+1,3)) * dk1;


    l1_now(j+1) = floor(si_b1(j+1,1) / sigma1);
    rho1_1(j+1) = si_b1(j+1,1) / sigma1  - l1_now(j+1);

    m1_now(j+1) = floor(si_b1(j+1,2) / sigma2);
    rho2_1(j+1) = si_b1(j+1,2) / sigma2  - m1_now(j+1);

    n1_now(j+1) = floor(si_b1(j+1,3) / sigma3);
    rho3_1(j+1) = si_b1(j+1,3) / sigma3  - n1_now(j+1);

    for i = 1 : j

        if l1_now(j+1) == l1_now(i)
            l1_num(j+1) = l1_num(i);
            break;
        end

        if i == j
            l_max = l_max + 1;
            l1_num(j+1) = l_max;
        end

    end

    for i = 1 : j

        if m1_now(j+1) == m1_now(i) 
            m1_num(j+1) = m1_num(i);
            break;
        end

        if i == j
            m_max = m_max + 1;
            m1_num(j+1) = m_max;
        end

    end

    for i = 1 : j

        if n1_now(j+1) == n1_now(i)
            n1_num(j+1) = n1_num(i);
            break;
        end

        if i == j
            n_max = n_max + 1;
            n1_num(j+1) = n_max;
        end

    end

    l = l1_num(j); %代入しやすくするため
    m = m1_num(j);
    n = n1_num(j);

    l2 = l1_num(j+1); %代入しやすくするため
    m2 = m1_num(j+1);
    n2 = n1_num(j+1);


    Ef1 = Ef1 + ( 0 - ( alpha(l,m,n) + rho1_1(j) * (alpha(l+1,m,n) - alpha(l,m,n)) + rho2_1(j) * (alpha(l,m+1,n) - alpha(l,m,n)) + rho3_1(j) * (alpha(l,m,n+1) - alpha(l,m,n))) ) ^ 2;

    Ef3 = Ef3 + ( y1(j) - (gamma(l,m,n) + rho1_1(j) * (gamma(l+1,m,n) - gamma(l,m,n)) + rho2_1(j) * (gamma(l,m+1,n) - gamma(l,m,n)) + rho3_1(j) * (gamma(l,m,n+1) - gamma(l,m,n))) ) ^ 2;

end

Ef1 = Ef1 + ( 0 - ( alpha(l2,m2,n2) + rho1_1(j+1) * (alpha(l2+1,m2,n2) - alpha(l2,m2,n2)) + rho2_1(j+1) * (alpha(l2,m2+1,n2) - alpha(l2,m2,n2)) + rho3_1(j+1) * (alpha(l2,m2,n2+1) - alpha(l2,m2,n2))) ) ^ 2;

Ef3 = Ef3 + ( y1(j+1) - (gamma(l2,m2,n2) + rho1_1(j+1) * (gamma(l2+1,m2,n2) - gamma(l2,m2,n2)) + rho2_1(j+1) * (gamma(l2,m2+1,n2) - gamma(l2,m2,n2)) + rho3_1(j+1) * (gamma(l2,m2,n2+1) - gamma(l2,m2,n2))) ) ^ 2;


%z1&z3軸の推定(k2軸上)----------------------------------------------------

dk2 = 1;   %時間刻み
K2fin = 14; %シミュレーション終了時間, length(k) = Kfin + 1
k2 = [0:dk2:K2fin];

u1_b2 = zeros(length(k2),1); %並進速度
u2_b2 = zeros(length(k2),1); %回転角速度


si_b2 = zeros(length(k2),3); %観測するセンサ変数 , s = (s1, s2, s3) = (x ,y, θ)
si_b2(1,:) = [1-1/sqrt(2) 3/sqrt(2) 3*pi/4];    %(s1, s2, s3)の初期値を設定

si_c2 = zeros(length(k2),3); %補正後のセンサ変数(zi,z3空間と等しい)、結果比較用
si_c2(1,:) = [1 1 pi/2];


f1_b2 = zeros(length(k2),1);   %写像f1:s→z1の推定
f3_b2 = zeros(length(k2),1);   %写像f3:s→z3の推定


l2_now = zeros(length(k2),1);  % l_now = 時刻kのl, 値が飛び飛び or 被る可能性あり, 0や負の値になる可能性あり
m2_now = zeros(length(k2),1);  % m_now = 時刻kのl, 値が飛び飛び or 被る可能性あり, 0や負の値になる可能性あり
n2_now = zeros(length(k2),1);  % n_now = 時刻kのl, 値が飛び飛び or 被る可能性あり, 0や負の値になる可能性あり

l2_num = ones(length(k2),1); %推定した格子点に前から順に1,2,3,...と番号をつけていく
m2_num = ones(length(k2),1);
n2_num = ones(length(k2),1);

rho1_2 = zeros(length(k2),1);
rho2_2 = zeros(length(k2),1);
rho3_2 = zeros(length(k2),1);


l2_now(1) = floor(si_b2(1,1) / sigma1);
rho1_2(1) = si_b2(1,1) / sigma1  - l2_now(1);

m2_now(1) = floor(si_b2(1,2) / sigma2);
rho2_2(1) = si_b2(1,2) / sigma2  - m2_now(1);

n2_now(1) = floor(si_b2(1,3) / sigma3);
rho3_2(1) = si_b2(1,3) / sigma3  - n2_now(1);


y2 = ones(length(k2),1);
y2(1) = 1;


for i = 1 : length(k2) - 1
        
    if l2_now(1) == l1_now(i)
        l2_num(1) = l1_num(i);
        break;
    end

    if i == length(k2) - 1
        l_max = l_max + 1;
        l2_num(1) = l_max;
    end

end

for i = 1 : length(k2) - 1
        
    if m2_now(1) == m1_now(i)
        m2_num(1) = m1_num(i);
        break;
    end

    if i == length(k2) - 1
        m_max = m_max + 1;
        m2_num(1) = m_max;
    end

end

for i = 1 : length(k2) - 1
        
    if n2_now(1) == n1_now(i)
        n2_num(1) = n1_num(i);
        break;
    end

    if i == length(k2) - 1
        n_max = n_max + 1;
        n2_num(1) = n_max;
    end

end

%---誤差関数の定義（k1軸上）----------------------------------------------------


for j = 1 : length(k2) - 1

    if mod(j,10) == 0
        u1_b2(j+1) = -1;
        u2_b2(j+1) = 0;
        y2(j+1) = y2(j) - 1;
    elseif mod(j,5) == 0
        u1_b2(j+1) = 1;
        u2_b2(j+1) = 0;
        y2(j+1) = y2(j) - 1;
    elseif j < 5
        u1_b2(j+1) = 0;
        u2_b2(j+1) = -pi/4;
        y2(j+1) = y2(j);
    elseif j > 5 && j < 10
        u1_b2(j+1) = 0;
        u2_b2(j+1) = pi/4;
        y2(j+1) = y2(j);
    else
        u1_b2(j+1) = 0;
        u2_b2(j+1) = -pi/4;
        y2(j+1) = y2(j);
    end

   
    si_b2(j+1,3) = si_b2(j,3) + u2_b2(j+1) * dk2;
    si_b2(j+1,1) = si_b2(j,1) + u1_b2(j+1) * cos(si_b2(j+1,3)) * dk2;
    si_b2(j+1,2) = si_b2(j,2) + u1_b2(j+1) * sin(si_b2(j+1,3)) * dk2;

    si_c2(j+1,3) = si_c2(j,3) + u2_b2(j+1) * dk2;
    si_c2(j+1,1) = si_c2(j,1) + u1_b2(j+1) * cos(si_c2(j+1,3)) * dk2;
    si_c2(j+1,2) = si_c2(j,2) + u1_b2(j+1) * sin(si_c2(j+1,3)) * dk2;


    l2_now(j+1) = floor(si_b2(j+1,1) / sigma1);
    rho1_2(j+1) = si_b2(j+1,1) / sigma1  - l2_now(j+1);

    m2_now(j+1) = floor(si_b2(j+1,2) / sigma2);
    rho2_2(j+1) = si_b2(j+1,2) / sigma2  - m2_now(j+1);

    n2_now(j+1) = floor(si_b2(j+1,3) / sigma3);
    rho3_2(j+1) = si_b2(j+1,3) / sigma3  - n2_now(j+1);


    for i = 1 : length(k2) - 1

        if i <= j && l2_now(j+1) == l2_now(i)
            l2_num(j+1) = l2_num(i);
            break;
        end
            
        if l2_now(j+1) == l1_now(i)
            l2_num(j+1) = l1_num(i);
            break;
        end

        if i == length(k2) - 1
            l_max = l_max + 1;
            l2_num(j+1) = l_max;
        end
    
    end

    for i = 1 : length(k2) - 1

        if i <= j && m2_now(j+1) == m2_now(i)
            m2_num(j+1) = m2_num(i);
            break;
        end

        if m2_now(j+1) == m1_now(i) 
            m2_num(j+1) = m1_num(i);
            break;
        end

        if i == length(k2) - 1
            m_max = m_max + 1;
            m2_num(j+1) = m_max;
        end

    end

    for i = 1 : length(k2) - 1

        if i <= j && n2_now(j+1) == n2_now(i)
            n2_num(j+1) = n2_num(i);
            break;
        end

        if n2_now(j+1) == n1_now(i)
            n2_num(j+1) = n1_num(i);
            break;
        end

        if i == length(k2) - 1
            n_max = n_max + 1;
            n2_num(j+1) = n_max;
        end

    end

    l = l2_num(j); %代入しやすくするため
    m = m2_num(j);
    n = n2_num(j);

    l2 = l2_num(j+1); %代入しやすくするため
    m2 = m2_num(j+1);
    n2 = n2_num(j+1);

    Ef1 = Ef1 + ( 1 - ( alpha(l,m,n) + rho1_2(j) * (alpha(l+1,m,n) - alpha(l,m,n)) + rho2_2(j) * (alpha(l,m+1,n) - alpha(l,m,n)) + rho3_2(j) * (alpha(l,m,n+1) - alpha(l,m,n))) ) ^ 2;

    Ef3 = Ef3 + ( y2(j) - (gamma(l,m,n) + rho1_2(j) * (gamma(l+1,m,n) - gamma(l,m,n)) + rho2_2(j) * (gamma(l,m+1,n) - gamma(l,m,n)) + rho3_2(j) * (gamma(l,m,n+1) - gamma(l,m,n))) ) ^ 2;

end

Ef1 = Ef1 + ( 1 - ( alpha(l2,m2,n2) + rho1_2(j+1) * (alpha(l2+1,m2,n2) - alpha(l2,m2,n2)) + rho2_2(j+1) * (alpha(l2,m2+1,n2) - alpha(l2,m2,n2)) + rho3_2(j+1) * (alpha(l2,m2,n2+1) - alpha(l2,m2,n2))) ) ^ 2;

Ef3 = Ef3 + ( y2(j+1) - (gamma(l2,m2,n2) + rho1_2(j+1) * (gamma(l2+1,m2,n2) - gamma(l2,m2,n2)) + rho2_2(j+1) * (gamma(l2,m2+1,n2) - gamma(l2,m2,n2)) + rho3_2(j+1) * (gamma(l2,m2,n2+1) - gamma(l2,m2,n2))) ) ^ 2;


%---正則化項の追加---------------

E4_f1 = 0;
E4_f3 = 0;

for a = 1 : l_max - 1
    for b = 1 : m_max - 1
        for c = 1 : n_max - 1

            E4_f1 = E4_f1 + (alpha(a+2,b,c) - 2 * alpha(a+1,b,c) - alpha(a,b,c)) ^ 2 + (alpha(a,b+2,c) - 2 * alpha(a,b+1,c) - alpha(a,b,c)) ^ 2 + (alpha(a,b,c+2) - 2 * alpha(a,b,c+1) - alpha(a,b,c)) ^ 2;
            E4_f3 = E4_f3 + (gamma(a+2,b,c) - 2 * gamma(a+1,b,c) + gamma(a,b,c)) ^ 2 + (gamma(a,b+2,c) - 2 * gamma(a,b+1,c) + gamma(a,b,c)) ^ 2 + (gamma(a,b,c+2) - 2 * gamma(a,b,c+1) + gamma(a,b,c)) ^ 2;

        end
    end
end

Ef1 = Ef1 + 0.2 * E4_f1;
Ef3 = Ef3 + 0.2 * E4_f3;



%---最急降下法によるパラメータの決定----------------------------

eta_f1 = 2.0 * 10 ^ (-1); %学習率
eta_f3 = 1.0 * 10 ^ (-1);

iteration = 300; %パラメータ更新回数（最大）

param_alpha = zeros(10,10,10,iteration+1);
param_gamma = zeros(10,10,10,iteration+1);



%---Eの設定---------------------------

Ef1_initial = double(subs(Ef1, [alpha(1:l_max+1,1:m_max+1,1:n_max+1)],[param_alpha(1:l_max+1,1:m_max+1,1:n_max+1,1)]));

Ef3_initial = double(subs(Ef3, [gamma(1:l_max+1,1:m_max+1,1:n_max+1)],[param_gamma(1:l_max+1,1:m_max+1,1:n_max+1,1)]));
                             

disp('Ef1 = ')
disp(Ef1_initial)
disp('Ef3 = ')
disp(Ef3_initial)
disp('--------------------')


Ef1_value = zeros(1,iteration);
Ef3_value = zeros(1,iteration);


for a = 1:l_max+1
    for b = 1:m_max+1
        for c = 1:n_max+1

            DEf1(a,b,c) = diff(Ef1,alpha(a,b,c));
            DEf3(a,b,c) = diff(Ef3,gamma(a,b,c));

            DEf1_1{a,b,c} = matlabFunction(DEf1(a,b,c), 'vars', {alpha(1:l_max+1,1:m_max+1,1:n_max+1)});
            DEf3_1{a,b,c} = matlabFunction(DEf3(a,b,c), 'vars', {gamma(1:l_max+1,1:m_max+1,1:n_max+1)});

            A = [a b c];
            disp(A)
            disp('------------')

        end
    end
end



for t = 1:iteration

    for a = 1:l_max+1
        for b = 1:m_max+1
            for c = 1:n_max+1

            DEf1_2 = DEf1_1{a,b,c}(param_alpha(1:l_max+1,1:m_max+1,1:n_max+1,t));
            DEf3_2 = DEf3_1{a,b,c}(param_gamma(1:l_max+1,1:m_max+1,1:n_max+1,t));

            param_alpha(a,b,c,t+1) = param_alpha(a,b,c,t) - eta_f1 * DEf1_2;
            param_gamma(a,b,c,t+1) = param_gamma(a,b,c,t) - eta_f3 * DEf3_2;

            end
        end
    end


    Ef1_value(t) = double(subs(Ef1, [alpha(1:l_max+1,1:m_max+1,1:n_max+1)],[param_alpha(1:l_max+1,1:m_max+1,1:n_max+1,t+1)]));
    
    Ef3_value(t) = double(subs(Ef3, [gamma(1:l_max+1,1:m_max+1,1:n_max+1)],[param_gamma(1:l_max+1,1:m_max+1,1:n_max+1,t+1)]));


    disp('t = ')
    disp(t)
    disp('Ef1 = ')
    disp(Ef1_value(t))
    disp('Ef3 = ')
    disp(Ef3_value(t))
 

    if t > 1

        if (Ef1_value(t) > Ef1_value(t-1))
            disp('Ef1_2が増加しました')
        end
        if (Ef3_value(t) > Ef3_value(t-1))
            disp('Ef3_2が増加しました')
        end
        if (Ef1_value(t) > Ef1_value(t-1)) || (Ef3_value(t) > Ef3_value(t-1))
            iteration = t;
            disp('iterationを強制終了します')
            break
        end
        if t == iteration
            iteration = t+1;
            disp('iterationを正常に終了することができました！')
            break
        end

    end

    disp('--------------------')

end

% 推定結果の取得--------------------------------------

for j = 1:length(k1)

    l = l1_num(j);
    m = m1_num(j);
    n = n1_num(j);

    f1_b1(j) = param_alpha(l,m,n,iteration) + rho1_1(j) * (param_alpha(l+1,m,n,iteration) - param_alpha(l,m,n,iteration))...
                                            + rho2_1(j) * (param_alpha(l,m+1,n,iteration) - param_alpha(l,m,n,iteration))...
                                            + rho3_1(j) * (param_alpha(l,m,n+1,iteration) - param_alpha(l,m,n,iteration));


    f3_b1(j) = param_gamma(l,m,n,iteration) + rho1_1(j) * (param_gamma(l+1,m,n,iteration) - param_gamma(l,m,n,iteration))...
                                            + rho2_1(j) * (param_gamma(l,m+1,n,iteration) - param_gamma(l,m,n,iteration))...
                                            + rho3_1(j) * (param_gamma(l,m,n+1,iteration) - param_gamma(l,m,n,iteration));

end


% 推定結果のplot--------------------------------------

figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c1(:,3), si_c1(:,1), '--m', si_c1(:,3), f1_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z1 = f1(s) = s1 の答え合わせ
xlabel("s3' = θ")
ylabel('z1 = f1(s)')
legend('真値：s1','推定値：z1 = f1(s)')

hold off;


figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c1(:,3), si_c1(:,2), '--m', si_c1(:,3), f3_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z3 = f3(s) = s2 の答え合わせ
xlabel("s3' = θ")
ylabel('z3 = f3(s)')
legend('真値：s2','推定値：z3 = f3(s)')

hold off;


% 推定結果の取得--------------------------------------

for j = 1:length(k2)

    l = l2_num(j);
    m = m2_num(j);
    n = n2_num(j);

    f1_b2(j) = param_alpha(l,m,n,iteration) + rho1_2(j) * (param_alpha(l+1,m,n,iteration) - param_alpha(l,m,n,iteration))...
                                            + rho2_2(j) * (param_alpha(l,m+1,n,iteration) - param_alpha(l,m,n,iteration))...
                                            + rho3_2(j) * (param_alpha(l,m,n+1,iteration) - param_alpha(l,m,n,iteration));


    f3_b2(j) = param_gamma(l,m,n,iteration) + rho1_2(j) * (param_gamma(l+1,m,n,iteration) - param_gamma(l,m,n,iteration))...
                                            + rho2_2(j) * (param_gamma(l,m+1,n,iteration) - param_gamma(l,m,n,iteration))...
                                            + rho3_2(j) * (param_gamma(l,m,n+1,iteration) - param_gamma(l,m,n,iteration));

end

toc

% 推定結果のplot--------------------------------------

figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c2(:,3), si_c2(:,1), '--m', si_c2(:,3), f1_b2(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z1 = f1(s) = s1 の答え合わせ
xlabel("s3' = θ")
ylabel('z1 = f1(s)')
legend('真値：s1','推定値：z1 = f1(s)')

hold off;


figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c2(:,3), si_c2(:,2), '--m', si_c2(:,3), f3_b2(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z3 = f3(s) = s2 の答え合わせ
xlabel("s3' = θ")
ylabel('z3 = f3(s)')
legend('真値：s2','推定値：z3 = f3(s)')

hold off;




%---f2,g,hの推定----------------------------------------------------

dk = 1;   %時間刻み
Kfin = 15; %シミュレーション終了時間
k = [0:dk:Kfin];

u1_b = ones(length(k),1) * 0.069;
u2_b = zeros(length(k),1);

si_b = zeros(length(k),3); %観測するセンサ変数 , s = (s1, s2, s3) = (x ,y, θ)
si_b(1,:) = [1-1/sqrt(2) 1+1/sqrt(2) 7*pi/12];   %(s1, s2, s3)=(x ,y, θ)の初期値を設定

si_c = zeros(length(k),3); %補正後のセンサ変数(zi,z3空間と等しい)、結果比較用
si_c(1,:) = [0 0 pi/3];

f2_b = zeros(length(k),1);   %写像f2:s→z2の推定
gmap_b = zeros(length(k),1);
hmap_b = zeros(length(k),1);

l_now = zeros(length(k),1);
m_now = zeros(length(k),1);
n_now = zeros(length(k),1);

l_num = ones(length(k),1);
m_num = ones(length(k),1);
n_num = ones(length(k),1);

rho1 = zeros(length(k),1);
rho2 = zeros(length(k),1);
rho3 = zeros(length(k),1);

beta = sym('beta',[10 10 10]);
delta = sym('delta',[10 10 10]);
epsilon = sym('epsilon',[10 10 10]);

Ef2 = 0;
Eg = 0;
Eh = 0;

l_now(1) = floor(si_b(1,1) / sigma1);
rho1(1) = si_b(1,1) / sigma1  - l_now(1);

m_now(1) = floor(si_b(1,2) / sigma2);
rho2(1) = si_b(1,2) / sigma2  - m_now(1);

n_now(1) = floor(si_b(1,3) / sigma3);
rho3(1) = si_b(1,3) / sigma3  - n_now(1);


for i = 1 : length(k1) - 1 %k1=k2より
        
    if l_now(1) == l1_now(i)
        l_num(1) = l1_num(i);
        break;
    end

    if l_now(1) == l2_now(i)
        l_num(1) = l2_num(i);
        break;
    end

    if i == length(k1) - 1
       break; %何かしら設定必要
    end

end

for i = 1 : length(k1) - 1
        
    if m_now(1) == m1_now(i)
        m_num(1) = m1_num(i);
        break;
    end

    if m_now(1) == m2_now(i)
        m_num(1) = m2_num(i);
        break;
    end

    if i == length(k1) - 1
       break; %何かしら設定必要
    end

end

for i = 1 : length(k1) - 1
        
    if n_now(1) == n1_now(i)
        n_num(1) = n1_num(i);
        break;
    end

    if n_now(1) == n2_now(i)
        n_num(1) = n2_num(i);
        break;
    end

    if i == length(k1) - 1
       break; %何かしら設定必要
    end

end


%---誤差関数の定義(f2, h)----------------------------------------------------

for j = 1:length(k) - 1

    si_b(j+1,3) = si_b(j,3) + u2_b(j) * dk;
    si_b(j+1,1) = si_b(j,1) + u1_b(j) * cos(si_b(j+1,3)) * dk;
    si_b(j+1,2) = si_b(j,2) + u1_b(j) * sin(si_b(j+1,3)) * dk;

    si_c(j+1,3) = si_c(j,3) + u2_b(j+1) * dk;
    si_c(j+1,1) = si_c(j,1) + u1_b(j+1) * cos(si_c(j+1,3)) * dk;
    si_c(j+1,2) = si_c(j,2) + u1_b(j+1) * sin(si_c(j+1,3)) * dk;


    l_now(j+1) = floor(si_b(j+1,1) / sigma1);
    rho1(j+1) = si_b(j+1,1) / sigma1  - l_now(j+1);

    m_now(j+1) = floor(si_b(j+1,2) / sigma2);
    rho2(j+1) = si_b(j+1,2) / sigma2  - m_now(j+1);

    n_now(j+1) = floor(si_b(j+1,3) / sigma3);
    rho3(j+1) = si_b(j+1,3) / sigma3  - n_now(j+1);


    for i = 1 : length(k) - 1

        if i <= j && l_now(j+1) == l_now(i)
            l_num(j+1) = l_num(i);
            break;
        end
            
        if i <= length(k1) && l_now(j+1) == l1_now(i)
            l_num(j+1) = l1_num(i);
            break;
        end

        if i <= length(k2) && l_now(j+1) == l2_now(i)
            l_num(j+1) = l2_num(i);
            break;
        end

        if i == length(k) - 1
            break; %何かしら設定必要
        end
    
    end

    for i = 1 : length(k) - 1

        if i <= j && m_now(j+1) == m_now(i)
            m_num(j+1) = m_num(i);
            break;
        end
            
        if i <= length(k1) && m_now(j+1) == m1_now(i)
            m_num(j+1) = m1_num(i);
            break;
        end

        if i <= length(k2) && m_now(j+1) == m2_now(i)
            m_num(j+1) = m2_num(i);
            break;
        end

        if i == length(k) - 1
            break; %何かしら設定必要
        end
    
    end

    for i = 1 : length(k) - 1

        if i <= j && n_now(j+1) == n_now(i)
            n_num(j+1) = n_num(i);
            break;
        end
            
        if i <= length(k1) && n_now(j+1) == n1_now(i)
            n_num(j+1) = n1_num(i);
            break;
        end

        if i <= length(k2) && n_now(j+1) == n2_now(i)
            n_num(j+1) = n2_num(i);
            break;
        end

        if i == length(k) - 1
            break; %何かしら設定必要
        end
    
    end


    l = l_num(j); %代入しやすくするため
    m = m_num(j);
    n = n_num(j);

    l2 = l_num(j+1); %代入しやすくするため
    m2 = m_num(j+1);
    n2 = n_num(j+1);

    Ef2 = Ef2 + ( beta(l,m,n) + rho1(j) * (beta(l+1,m,n) - beta(l,m,n)) + rho2(j) * (beta(l,m+1,n) - beta(l,m,n)) + rho3(j) * (beta(l,m,n+1) - beta(l,m,n))...
                 - ( (param_gamma(l2,m2,n2,iteration) + rho1(j+1) * (param_gamma(l2+1,m2,n2,iteration) - param_gamma(l2,m2,n2,iteration)) + rho2(j+1) * (param_gamma(l2,m2+1,n2,iteration) - param_gamma(l2,m2,n2,iteration)) + rho3(j+1) * (param_gamma(l2,m2,n2+1,iteration) - param_gamma(l2,m2,n2,iteration)))...
                 - (param_gamma(l,m,n,iteration) + rho1(j) * (param_gamma(l+1,m,n,iteration) - param_gamma(l,m,n,iteration)) + rho2(j) * (param_gamma(l,m+1,n,iteration) - param_gamma(l,m,n,iteration)) + rho3(j) * (param_gamma(l,m,n+1,iteration) - param_gamma(l,m,n,iteration))) )...
                 / ( (param_alpha(l2,m2,n2,iteration) + rho1(j+1) * (param_alpha(l2+1,m2,n2,iteration) - param_alpha(l2,m2,n2,iteration)) + rho2(j+1) * (param_alpha(l2,m2+1,n2,iteration) - param_alpha(l2,m2,n2,iteration)) + rho3(j+1) * (param_alpha(l2,m2,n2+1,iteration) - param_alpha(l2,m2,n2,iteration)))...
                 - (param_alpha(l,m,n,iteration) + rho1(j) * (param_alpha(l+1,m,n,iteration) - param_alpha(l,m,n,iteration)) + rho2(j) * (param_alpha(l,m+1,n,iteration) - param_alpha(l,m,n,iteration)) + rho3(j) * (param_alpha(l,m,n+1,iteration) - param_alpha(l,m,n,iteration))) ) ) ^ 2;

    % E2 = E2 + ( u2_b(j) / u1_b(j) * ( delta(l,m,n) + rho1(j) * (delta(l+1,m,n) - delta(l,m,n)) + rho2(j) * (delta(l,m+1,n) - delta(l,m,n)) + rho3(j) * (delta(l,m,n+1) - delta(l,m,n)) )...
    %              - ( (beta(l2,m2,n2) + rho1(j+1) * (beta(l2+1,m2,n2) - beta(l2,m2,n2)) + rho2(j+1) * (beta(l2,m2+1,n2) - beta(l2,m2,n2)) + rho3(j+1) * (beta(l2,m2,n2+1) - beta(l2,m2,n2)))...
    %              - (beta(l,m,n) + rho1(j) * (beta(l+1,m,n) - beta(l,m,n)) + rho2(j) * (beta(l,m+1,n) - beta(l,m,n)) + rho3(j) * (beta(l,m,n+1) - beta(l,m,n))) )...
    %              / ( (alpha(l2,m2,n2) + rho1(j+1) * (alpha(l2+1,m2,n2) - alpha(l2,m2,n2)) + rho2(j+1) * (alpha(l2,m2+1,n2) - alpha(l2,m2,n2)) + rho3(j+1) * (alpha(l2,m2,n2+1) - alpha(l2,m2,n2)))...
    %              - (alpha(l,m,n) + rho1(j) * (alpha(l+1,m,n) - alpha(l,m,n)) + rho2(j) * (alpha(l,m+1,n) - alpha(l,m,n)) + rho3(j) * (alpha(l,m,n+1) - alpha(l,m,n))) ) ) ^ 2;
    
    Eh = Eh + ( u1_b(j) * ( epsilon(l,m,n) + rho1(j) * (epsilon(l+1,m,n) - epsilon(l,m,n)) + rho2(j) * (epsilon(l,m+1,n) - epsilon(l,m,n)) + rho3(j) * (epsilon(l,m,n+1) - epsilon(l,m,n)) )...
                 - ( (param_alpha(l2,m2,n2,iteration) + rho1(j+1) * (param_alpha(l2+1,m2,n2,iteration) - param_alpha(l2,m2,n2,iteration)) + rho2(j+1) * (param_alpha(l2,m2+1,n2,iteration) - param_alpha(l2,m2,n2,iteration)) + rho3(j+1) * (param_alpha(l2,m2,n2+1,iteration) - param_alpha(l2,m2,n2,iteration)))...
                 - (param_alpha(l,m,n,iteration) + rho1(j) * (param_alpha(l+1,m,n,iteration) - param_alpha(l,m,n,iteration)) + rho2(j) * (param_alpha(l,m+1,n,iteration) - param_alpha(l,m,n,iteration)) + rho3(j) * (param_alpha(l,m,n+1,iteration) - param_alpha(l,m,n,iteration))) )...
                 / dk ) ^ 2;

end


%---最急降下法によるパラメータの決定----------------------------

eta_f2 = 1.0 * 10 ^ (-1); %学習率
eta_h = 50 * 10 ^ (-1);

iteration = 300; %パラメータ更新回数（最大）

param_beta = zeros(10,10,10,iteration+1);
param_epsilon = zeros(10,10,10,iteration+1);


%---Eの設定---------------------------

Ef2_initial = double(subs(Ef2, [beta(1:l_max+1,1:m_max+1,1:n_max+1)],[param_beta(1:l_max+1,1:m_max+1,1:n_max+1,1)]));

Eh_initial = double(subs(Eh, [epsilon(1:l_max+1,1:m_max+1,1:n_max+1)],[param_epsilon(1:l_max+1,1:m_max+1,1:n_max+1,1)]));
                             

disp('Ef2 = ')
disp(Ef2_initial)
disp('Eh = ')
disp(Eh_initial)
disp('--------------------')


Ef2_value = zeros(1,iteration);
Eh_value = zeros(1,iteration);


for a = 1:l_max+1
    for b = 1:m_max+1
        for c = 1:n_max+1

            DEf2(a,b,c) = diff(Ef2,beta(a,b,c));
            DEh(a,b,c) = diff(Eh,epsilon(a,b,c));

            DEf2_1{a,b,c} = matlabFunction(DEf2(a,b,c), 'vars', {beta(1:l_max+1,1:m_max+1,1:n_max+1)});
            DEh_1{a,b,c} = matlabFunction(DEh(a,b,c), 'vars', {epsilon(1:l_max+1,1:m_max+1,1:n_max+1)});

            A = [a b c];
            disp(A)
            disp('------------')

        end
    end
end



for t = 1:iteration

    for a = 1:l_max+1
        for b = 1:m_max+1
            for c = 1:n_max+1

            DEf2_2 = DEf2_1{a,b,c}(param_beta(1:l_max+1,1:m_max+1,1:n_max+1,t));
            DEh_2 = DEh_1{a,b,c}(param_epsilon(1:l_max+1,1:m_max+1,1:n_max+1,t));

            param_beta(a,b,c,t+1) = param_beta(a,b,c,t) - eta_f2 * DEf2_2;
            param_epsilon(a,b,c,t+1) = param_epsilon(a,b,c,t) - eta_h * DEh_2;

            end
        end
    end


    Ef2_value(t) = double(subs(Ef2, [beta(1:l_max+1,1:m_max+1,1:n_max+1)],[param_beta(1:l_max+1,1:m_max+1,1:n_max+1,t+1)]));
    
    Eh_value(t) = double(subs(Eh, [epsilon(1:l_max+1,1:m_max+1,1:n_max+1)],[param_epsilon(1:l_max+1,1:m_max+1,1:n_max+1,t+1)]));


    disp('t = ')
    disp(t)
    disp('Ef2 = ')
    disp(Ef2_value(t))
    disp('Eh = ')
    disp(Eh_value(t))
 

    if t > 1

        if (Ef2_value(t) > Ef2_value(t-1))
            disp('Ef2が増加しました')
        end
        if (Eh_value(t) > Eh_value(t-1))
            disp('Ehが増加しました')
        end
        if (Ef2_value(t) > Ef2_value(t-1)) || (Eh_value(t) > Eh_value(t-1))
            iteration = t;
            disp('iterationを強制終了します')
            break
        end
        if t == iteration
            iteration = t+1;
            disp('iterationを正常に終了することができました！')
            break
        end

    end

    disp('--------------------')

end


% 推定結果の取得--------------------------------------

for j = 1:length(k)

    l = l_num(j);
    m = m_num(j);
    n = n_num(j);

    f1_b(j) = param_alpha(l,m,n,iteration) + rho1(j) * (param_alpha(l+1,m,n,iteration) - param_alpha(l,m,n,iteration))...
                                           + rho2(j) * (param_alpha(l,m+1,n,iteration) - param_alpha(l,m,n,iteration))...
                                           + rho3(j) * (param_alpha(l,m,n+1,iteration) - param_alpha(l,m,n,iteration));

    f3_b(j) = param_beta(l,m,n,iteration) + rho1(j) * (param_gamma(l+1,m,n,iteration) - param_gamma(l,m,n,iteration))...
                                          + rho2(j) * (param_gamma(l,m+1,n,iteration) - param_gamma(l,m,n,iteration))...
                                          + rho3(j) * (param_gamma(l,m,n+1,iteration) - param_gamma(l,m,n,iteration));


    f2_b(j) = param_beta(l,m,n,iteration) + rho1(j) * (param_beta(l+1,m,n,iteration) - param_beta(l,m,n,iteration))...
                                          + rho2(j) * (param_beta(l,m+1,n,iteration) - param_beta(l,m,n,iteration))...
                                          + rho3(j) * (param_beta(l,m,n+1,iteration) - param_beta(l,m,n,iteration));


    hmap_b(j) = param_epsilon(l,m,n,iteration) + rho1(j) * (param_epsilon(l+1,m,n,iteration) - param_epsilon(l,m,n,iteration))...
                                            + rho2(j) * (param_epsilon(l,m+1,n,iteration) - param_epsilon(l,m,n,iteration))...
                                            + rho3(j) * (param_epsilon(l,m,n+1,iteration) - param_epsilon(l,m,n,iteration));

end


% 推定結果のplot--------------------------------------

figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c(:,3), si_c(:,1), '--m', si_c(:,3), f1_b(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z1 = f1(s) = s1 の答え合わせ
xlabel("s3' = θ")
ylabel('z1 = f1(s)')
legend('真値：s1','推定値：z1 = f1(s)')

hold off;


figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c(:,3), si_c(:,2), '--m', si_c(:,3), f3_b(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z3 = f3(s) = s2 の答え合わせ
xlabel("s3' = θ")
ylabel('z3 = f3(s)')
legend('真値：s2','推定値：z3 = f3(s)')

hold off;




% 推定結果のplot--------------------------------------

figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c(:,3), tan(si_c(:,3)), '--m', si_c(:,3), f2_b(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z1 = f1(s) = s1 の答え合わせ
xlabel("s3' = θ")
ylabel('z2 = f2(s)')
legend('真値：tan(s3)','推定値：z2 = f2(s)')

hold off;


figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c(:,3), cos(si_c(:,3)), '--m', si_c(:,3), hmap_b(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z3 = f3(s) = s2 の答え合わせ
xlabel("s3' = θ")
ylabel('h(s)')
legend('真値：cos(s3)','推定値：h(s)')

hold off;
