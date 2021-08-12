clear;
close all;

tic

%z1&z3軸の推定(k1軸上)----------------------------------------------------

dk1 = 1;   %時間刻み
K1fin = 9; %シミュレーション終了時間, length(k) = Kfin + 1
k1 = [0:dk1:K1fin];

u1_b1 = zeros(length(k1),1); %並進速度
u2_b1 = zeros(length(k1),1); %回転角速度


si_b1 = zeros(length(k1),3); %観測するセンサ変数 , s = (s1, s2, s3) = (x ,y, θ)
si_b1(1,:) = [1 1 -pi/4];    %(s1, s2, s3)の初期値を設定

si_c1 = zeros(length(k1),3); %補正後のセンサ変数 , s' = (s1', s2', s3') = (s1　-　s1_initial,s2 - s2_initial,s3 - s3_initial)
si_c1(1,:) = [0 0 -pi/2];    %(s1, s2, s3)の初期値を設定


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

alpha = sym('alpha',[50 50 50]); %l,m,nの順
beta = sym('beta',[50 50 50]);
gamma = sym('gamma',[50 50 50]);
delta = sym('delta',[50 50 50]);
epsilon = sym('epsilon',[50 50 50]);


sigma1 = 1; %s1のスケーリング定数
sigma2 = 1; %s1のスケーリング定数
sigma3 = pi/4; %s1のスケーリング定数

Ef1_1 = 0;
Ef3_1 = 0;
E4_f1_1 = 0;
E4_f3_1 = 0;

l1_now(1) = floor(si_b1(1,1) / sigma1);
rho1_1(1) = si_b1(1,1) / sigma1  - l1_now(1);

m1_now(1) = floor(si_b1(1,2) / sigma2);
rho2_1(1) = si_b1(1,2) / sigma2  - m1_now(1);

n1_now(1) = floor(si_b1(1,3) / sigma3);
rho3_1(1) = si_b1(1,3) / sigma3  - n1_now(1);


%---誤差関数の定義（k1軸上）----------------------------------------------------


for j = 1 : length(k1) - 1

    if mod(j,5) == 0
        u1_b1(j+1) = 1;
        u2_b1(j+1) = 0;
    elseif j < 5
        u1_b1(j+1) = 0;
        u2_b1(j+1) = pi/4;
    else
        u1_b1(j+1) = 0;
        u2_b1(j+1) = -pi/4;
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


    if l1_now(j+1) == l1_now(j)
        l1_num(j+1) = l1_num(j);
    else
        l1_num(j+1) = l1_num(j) + 1; % 対応する格子点に番号をつけていく
    end

    if m1_now(j+1) == m1_now(j) 
        m1_num(j+1) = m1_num(j);
    else
        m1_num(j+1) = m1_num(j) + 1; % 対応する格子点に番号をつけていく
    end

    if n1_now(j+1) == n1_now(j)
        n1_num(j+1) = n1_num(j);
    else
        n1_num(j+1) = n1_num(j) + 1; % 対応する格子点に番号をつけていく
    end

    l = l1_num(j); %代入しやすくするため
    m = m1_num(j);
    n = n1_num(j);

    l2 = l1_num(j+1); %代入しやすくするため
    m2 = m1_num(j+1);
    n2 = n1_num(j+1);


    Ef1_1 = Ef1_1 + ( 0 - ( (alpha(l2,m2,n2) + rho1_1(j+1) * (alpha(l2+1,m2,n2) - alpha(l2,m2,n2)) + rho2_1(j+1) * (alpha(l2,m2+1,n2) - alpha(l2,m2,n2)) + rho3_1(j+1) * (alpha(l2,m2,n2+1) - alpha(l2,m2,n2)))...
                        - (alpha(l,m,n) + rho1_1(j) * (alpha(l+1,m,n) - alpha(l,m,n)) + rho2_1(j) * (alpha(l,m+1,n) - alpha(l,m,n)) + rho3_1(j) * (alpha(l,m,n+1) - alpha(l,m,n))) )...
                        / dk1 ) ^ 2;

    Ef3_1 = Ef3_1 + ( u1_b1(j+1) - ( (gamma(l2,m2,n2) + rho1_1(j+1) * (gamma(l2+1,m2,n2) - gamma(l2,m2,n2)) + rho2_1(j+1) * (gamma(l2,m2+1,n2) - gamma(l2,m2,n2)) + rho3_1(j+1) * (gamma(l2,m2,n2+1) - gamma(l2,m2,n2)))...
                        - ( gamma(l,m,n) + rho1_1(j) * (gamma(l+1,m,n) - gamma(l,m,n)) + rho2_1(j) * (gamma(l,m+1,n) - gamma(l,m,n)) + rho3_1(j) * (gamma(l,m,n+1) - gamma(l,m,n))) )...
                        / dk1 ) ^ 2;


end



%---最急降下法によるパラメータの決定----------------------------

eta_f1 = 1.0 * 10 ^ (-1); %学習率
eta_f3 = 2.5 * 10 ^ (-1);

iteration = 300; %パラメータ更新回数（最大）

param_alpha = zeros(l2+1,m2+1,n2+1,iteration+1);
param_gamma = zeros(l2+1,m2+1,n2+1,iteration+1);


%---Eの設定-----------------------------

Ef1_1_initial = double(subs(Ef1_1, [alpha(1:l2+1,1:m2+1,1:n2+1)],[param_alpha(:,:,:,1)]));

Ef3_1_initial = double(subs(Ef3_1, [gamma(1:l2+1,1:m2+1,1:n2+1)],[param_gamma(:,:,:,1)]));
                             

disp('Ef1_1 = ')
disp(Ef1_1_initial)
disp('Ef3_1 = ')
disp(Ef3_1_initial)
disp('--------------------')


Ef1_1_value = zeros(1,iteration);
Ef3_1_value = zeros(1,iteration);


for a = 1:l2+1
    for b = 1:m2+1
        for c = 1:n2+1

            DEf1(a,b,c) = diff(Ef1_1,alpha(a,b,c));
            DEf3(a,b,c) = diff(Ef3_1,gamma(a,b,c));

            DEf1_1{a,b,c} = matlabFunction(DEf1(a,b,c), 'vars', {alpha(1:l2+1,1:m2+1,1:n2+1)});
            DEf3_1{a,b,c} = matlabFunction(DEf3(a,b,c), 'vars', {gamma(1:l2+1,1:m2+1,1:n2+1)});

            A = [a b c];
            disp(A)
            disp('------------')

        end
    end
end



for t = 1:iteration

    for a = 1:l2+1
        for b = 1:m2+1
            for c = 1:n2+1

            DEf1_2 = DEf1_1{a,b,c}(param_alpha(:,:,:,t));
            DEf3_2 = DEf3_1{a,b,c}(param_gamma(:,:,:,t));

            param_alpha(a,b,c,t+1) = param_alpha(a,b,c,t) - eta_f1 * DEf1_2;
            param_gamma(a,b,c,t+1) = param_gamma(a,b,c,t) - eta_f3 * DEf3_2;

            end
        end
    end

    param_alpha(1,1,1,t+1) = 0;
    param_gamma(1,1,1,t+1) = 0;


    Ef1_1_value(t) = double(subs(Ef1_1, [alpha(1:l2+1,1:m2+1,1:n2+1)],[param_alpha(:,:,:,t+1)]));
    
    Ef3_1_value(t) = double(subs(Ef3_1, [gamma(1:l2+1,1:m2+1,1:n2+1)],[param_gamma(:,:,:,t+1)]));


    disp('t = ')
    disp(t)
    disp('Ef1_1 = ')
    disp(Ef1_1_value(t))
    disp('Ef3 = ')
    disp(Ef3_1_value(t))
 

    if t > 1

        if (Ef1_1_value(t) > Ef1_1_value(t-1))
            disp('Ef1_1が増加しました')
        end
        if (Ef3_1_value(t) > Ef3_1_value(t-1))
            disp('Ef3_1が増加しました')
        end
        if (Ef1_1_value(t) > Ef1_1_value(t-1)) || (Ef3_1_value(t) > Ef3_1_value(t-1))
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

toc

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

