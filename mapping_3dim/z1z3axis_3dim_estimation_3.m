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

si_c1 = zeros(length(k1),3); %補正後のセンサ変数 , s' = (s1', s2', s3') = (s1　-　s1_initial,s2 - s2_initial,s3 - s3_initial)
si_c1(1,:) = [0 -1 -pi/2];    %(s1, s2, s3)の初期値を設定


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
% beta = sym('beta',[50 50 50]);
gamma = sym('gamma',[10 10 10]);
% delta = sym('delta',[50 50 50]);
% epsilon = sym('epsilon',[50 50 50]);


sigma1 = 1; %s1のスケーリング定数
sigma2 = 1; %s1のスケーリング定数
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
si_b2(1,:) = [1 1+sqrt(2) 3*pi/4];    %(s1, s2, s3)の初期値を設定

si_c2 = zeros(length(k2),3); %補正後のセンサ変数 , s' = (s1', s2', s3') = (s1　-　s1_initial,s2 - s2_initial,s3 - s3_initial)
si_c2(1,:) = [1 1 pi/2];    %(s1, s2, s3)の初期値を設定


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



%---最急降下法によるパラメータの決定----------------------------

eta_f1 = 2.5 * 10 ^ (-1); %学習率
eta_f3 = 2.5 * 10 ^ (-1);

iteration = 50; %パラメータ更新回数（最大）

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



