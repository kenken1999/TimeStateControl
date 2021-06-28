clear;
close all;

%all(f1,f2,f3,g,h)_estimation----------------------------------------------------

dk = 0.02;   %時間刻み
Kfin = 0.10; %シミュレーション終了時間
k = [0:dk:Kfin];

u1_b = ones(1,length(k)) * 5;
u2_b = ones(1,length(k)) * 5;

si_b = zeros(length(k),3); %観測するセンサ変数 , s = (s1, s2, s3) = (x ,y, θ)
si_b(1,:) = [0 0 -pi/2];   %(s1, s2, s3)=(x ,y, θ)の初期値を設定

f1_b = zeros(1,length(k));   %写像f1:s→z1の推定
f2_b = zeros(1,length(k));   %写像f2:s→z2の推定
f3_b = zeros(1,length(k));   %写像f3:s→z3の推定
gmap_b = zeros(1,length(k)); %dz2/dz1 = μ2 = u2/u1 * g(s)の推定, g(s) = 1/cos(s3)^3
hmap_b = zeros(1,length(k)); %dz1/dt = μ1 = u1*h(s)の推定, h(s) = cos(s3)

l_now = zeros(1,length(k));  % l_now = 時刻kのl, 値が飛び飛び or 被る可能性あり
m_now = zeros(1,length(k));  % m_now = 時刻kのl, 値が飛び飛び or 被る可能性あり
n_now = zeros(1,length(k));  % n_now = 時刻kのl, 値が飛び飛び or 被る可能性あり

l_num = ones(1,length(k)); %推定した格子点に前から順に番号をつけていく
m_num = ones(1,length(k));
n_num = ones(1,length(k));

alpha = sym('alpha',[50 50 50]); %l,m,nの順
beta = sym('beta',[50 50 50]);
gamma = sym('gamma',[50 50 50]);
delta = sym('delta',[50 50 50]);
epsilon = sym('epsilon',[50 50 50]);

sigma1 = 0.1; %s1のスケーリング定数
sigma2 = 0.1; %s2のスケーリング定数
sigma3 = 0.1; %s3のスケーリング定数

zeta1 = 0.1; %E1の調整係数
zeta2 = 0.1; %E2の調整係数
zeta3 = 0.1; %E3の調整係数

E = 0;

for j = 1:length(k) - 1

    si_b(j+1,3) = si_b(j,3) + u2_b(j) * dk;
    si_b(j+1,1) = si_b(j,1) + u1_b(j) * cos(si_b(j+1,3)) * dk;
    si_b(j+1,2) = si_b(j,2) + u1_b(j) * sin(si_b(j+1,3)) * dk;

    l_now(j+1) = floor(si_b(j+1,1) / sigma1);
    rho1 = si_b(j+1,1) / sigma1  - l_now(j+1);

    m_now(j+1) = floor(si_b(j+1,2) / sigma2);
    rho2 = si_b(j+1,2) / sigma2  - m_now(j+1);

    n_now(j+1) = floor(si_b(j+1,3) / sigma3);
    rho3 = si_b(j+1,3) / sigma3  - n_now(j+1);

    if l_now(j+1) > l_now(j)
        l_num(j+1) = l_num(j) + 1; % 対応する格子点に番号をつけていく
    else
        l_num(j+1) = l_num(j);
    end

    if m_now(j+1) > m_now(j)
        m_num(j+1) = m_num(j) + 1; % 対応する格子点に番号をつけていく
    else
        m_num(j+1) = m_num(j);
    end

    if n_now(j+1) > n_now(j)
        n_num(j+1) = n_num(j) + 1; % 対応する格子点に番号をつけていく
    else
        n_num(j+1) = n_num(j);
    end

    l = l_num(j); %代入しやすくするため
    m = m_num(j);
    n = n_num(j);

    l2 = l_num(j+1); %代入しやすくするため
    m2 = m_num(j+1);
    n2 = n_num(j+1);

    % E1 = 3 * beta(l,m,n) + rho1 * (beta(l+1,m,n) - beta(l,m,n)) + rho2 * (beta(l,m+1,n) - beta(l,m,n)) + rho3 * (beta(l,m,n+1) - beta(l,m,n)) - ( (3 * gamma(l2,m2,n2) + rho1 * (gamma(l2+1,m2,n2) - gamma(l2,m2,n2)) + rho2 * (gamma(l2,m2+1,n2) - gamma(l2,m2,n2)) + rho3 * (gamma(l2,m2,n2+1) - gamma(l2,m2,n2))) - (3 * gamma(l,m,n) + rho1 * (gamma(l+1,m,n) - gamma(l,m,n)) + rho2 * (gamma(l,m+1,n) - gamma(l,m,n)) + rho3 * (gamma(l,m,n+1) - gamma(l,m,n))) )  / ( (3 * alpha(l2,m2,n2) + rho1 * (alpha(l2+1,m2,n2) - alpha(l2,m2,n2)) + rho2 * (alpha(l2,m2+1,n2) - alpha(l2,m2,n2)) + rho3 * (alpha(l2,m2,n2+1) - alpha(l2,m2,n2))) - (3 * alpha(l,m,n) + rho1 * (alpha(l+1,m,n) - alpha(l,m,n)) + rho2 * (alpha(l,m+1,n) - alpha(l,m,n)) + rho3 * (alpha(l,m,n+1) - alpha(l,m,n))) );

    E1 = 3 * beta(l,m,n) + rho1 * (beta(l+1,m,n) - beta(l,m,n)) + rho2 * (beta(l,m+1,n) - beta(l,m,n)) + rho3 * (beta(l,m,n+1) - beta(l,m,n))...
             - ( (3 * gamma(l2,m2,n2) + rho1 * (gamma(l2+1,m2,n2) - gamma(l2,m2,n2)) + rho2 * (gamma(l2,m2+1,n2) - gamma(l2,m2,n2)) + rho3 * (gamma(l2,m2,n2+1) - gamma(l2,m2,n2)))...
             - (3 * gamma(l,m,n) + rho1 * (gamma(l+1,m,n) - gamma(l,m,n)) + rho2 * (gamma(l,m+1,n) - gamma(l,m,n)) + rho3 * (gamma(l,m,n+1) - gamma(l,m,n))) )...
         / ( (3 * alpha(l2,m2,n2) + rho1 * (alpha(l2+1,m2,n2) - alpha(l2,m2,n2)) + rho2 * (alpha(l2,m2+1,n2) - alpha(l2,m2,n2)) + rho3 * (alpha(l2,m2,n2+1) - alpha(l2,m2,n2)))...
             - (3 * alpha(l,m,n) + rho1 * (alpha(l+1,m,n) - alpha(l,m,n)) + rho2 * (alpha(l,m+1,n) - alpha(l,m,n)) + rho3 * (alpha(l,m,n+1) - alpha(l,m,n))) );

    E2 = u2_b(j) / u1_b(j) * ( 3 * delta(l,m,n) + rho1 * (delta(l+1,m,n) - delta(l,m,n)) + rho2 * (delta(l,m+1,n) - delta(l,m,n)) + rho3 * (delta(l,m,n+1) - delta(l,m,n)) )...
    - ( (3 * beta(l2,m2,n2) + rho1 * (beta(l2+1,m2,n2) - beta(l2,m2,n2)) + rho2 * (beta(l2,m2+1,n2) - beta(l2,m2,n2)) + rho3 * (beta(l2,m2,n2+1) - beta(l2,m2,n2)))...
        - (3 * beta(l,m,n) + rho1 * (beta(l+1,m,n) - beta(l,m,n)) + rho2 * (beta(l,m+1,n) - beta(l,m,n)) + rho3 * (beta(l,m,n+1) - beta(l,m,n))) )...
    / ( (3 * alpha(l2,m2,n2) + rho1 * (alpha(l2+1,m2,n2) - alpha(l2,m2,n2)) + rho2 * (alpha(l2,m2+1,n2) - alpha(l2,m2,n2)) + rho3 * (alpha(l2,m2,n2+1) - alpha(l2,m2,n2)))...
        - (3 * alpha(l,m,n) + rho1 * (alpha(l+1,m,n) - alpha(l,m,n)) + rho2 * (alpha(l,m+1,n) - alpha(l,m,n)) + rho3 * (alpha(l,m,n+1) - alpha(l,m,n))) );
    
    E3 = u1_b(j) * ( 3 * epsilon(l,m,n) + rho1 * (epsilon(l+1,m,n) - epsilon(l,m,n)) + rho2 * (epsilon(l,m+1,n) - epsilon(l,m,n)) + rho3 * (epsilon(l,m,n+1) - epsilon(l,m,n)) )...
    - ( (3 * alpha(l2,m2,n2) + rho1 * (alpha(l2+1,m2,n2) - alpha(l2,m2,n2)) + rho2 * (alpha(l2,m2+1,n2) - alpha(l2,m2,n2)) + rho3 * (alpha(l2,m2,n2+1) - alpha(l2,m2,n2)))...
        - (3 * alpha(l,m,n) + rho1 * (alpha(l+1,m,n) - alpha(l,m,n)) + rho2 * (alpha(l,m+1,n) - alpha(l,m,n)) + rho3 * (alpha(l,m,n+1) - alpha(l,m,n))) ) / dk;


    E = E + zeta1 * E1 ^ 2 + zeta2 * E2 ^ 2 + zeta3 * E3 ^ 2;
    % E = E + zeta1 * E1;


end

disp('l_now = ')
disp(l_now)

disp('m_now = ')
disp(m_now)

disp('n_now = ')
disp(n_now)

disp('l = ')
disp(l)

disp('m = ')
disp(m)

disp('n = ')
disp(n)


%---正則化項の追加---------------------------------
% for i = 2:p

%     Ef = Ef + (alpha(i+1) - 2 * alpha(i) + alpha(i-1)) ^ 2;

% end


%---写像 f1,f2,f3,g,h の推定----------------------------

eta_f1 = 0.05; %学習率
eta_f2 = 0.05;
eta_f3 = 0.05;
eta_g = 0.5;
eta_h = 0.01;

iteration = 3; %パラメータ更新回数（最大）

param_alpha = zeros(l+1,m+1,n+1,iteration);
param_beta = zeros(l+1,m+1,n+1,iteration);
param_gamma = zeros(l+1,m+1,n+1,iteration);
param_delta = zeros(l+1,m+1,n+1,iteration);
param_epsilon = zeros(l+1,m+1,n+1,iteration);

E_value = zeros(1,iteration);

disp('ok')

for t = 1:iteration

    for a = 1:l+1
        for b = 1:m+1
            for c = 1:n+1

                DEf1 = diff(E,alpha(a,b,c));
                DEf2 = diff(E,beta(a,b,c));
                DEf3 = diff(E,gamma(a,b,c));
                DEg = diff(E,delta(a,b,c));
                DEh = diff(E,epsilon(a,b,c));

                disp('good1')

                DEf1_2 = subs(DEf1, [alpha(1:l+1,1:m+1,1:n+1),beta(1:l+1,1:m+1,1:n+1),gamma(1:l+1,1:m+1,1:n+1),delta(1:l+1,1:m+1,1:n+1),epsilon(1:l+1,1:m+1,1:n+1)],...
                                    [param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t)]);

                disp('good2')
            

                DEf2_2 = subs(DEf2, [alpha(1:l+1,1:m+1,1:n+1),beta(1:l+1,1:m+1,1:n+1),gamma(1:l+1,1:m+1,1:n+1),delta(1:l+1,1:m+1,1:n+1),epsilon(1:l+1,1:m+1,1:n+1)],...
                                    [param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t)]);

                disp('good3')

                DEf3_2 = subs(DEf3, [alpha(1:l+1,1:m+1,1:n+1),beta(1:l+1,1:m+1,1:n+1),gamma(1:l+1,1:m+1,1:n+1),delta(1:l+1,1:m+1,1:n+1),epsilon(1:l+1,1:m+1,1:n+1)],...
                                    [param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t)]);

                DEg_2 = subs(DEg, [alpha(1:l+1,1:m+1,1:n+1),beta(1:l+1,1:m+1,1:n+1),gamma(1:l+1,1:m+1,1:n+1),delta(1:l+1,1:m+1,1:n+1),epsilon(1:l+1,1:m+1,1:n+1)],...
                                    [param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t)]);


                DEh_2 = subs(DEh, [alpha(1:l+1,1:m+1,1:n+1),beta(1:l+1,1:m+1,1:n+1),gamma(1:l+1,1:m+1,1:n+1),delta(1:l+1,1:m+1,1:n+1),epsilon(1:l+1,1:m+1,1:n+1)],...
                                    [param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t)]);

                param_alpha(a,b,c,t+1) = param_alpha(a,b,c,t) - eta_f1 * double(DEf1_2);
                param_beta(a,b,c,t+1) = param_beta(a,b,c,t) - eta_f2 * double(DEf2_2);
                param_gamma(a,b,c,t+1) = param_gamma(a,b,c,t) - eta_f3 * double(DEf3_2);
                param_delta(a,b,c,t+1) = param_delta(a,b,c,t) - eta_g * double(DEg_2);
                param_epsilon(a,b,c,t+1) = param_epsilon(a,b,c,t) - eta_h * double(DEh_2);

                disp('good4')

            end
        end
    end

    E_value(t) = double(subs(E, [alpha(1:l+1,1:m+1,1:n+1),beta(1:l+1,1:m+1,1:n+1),gamma(1:l+1,1:m+1,1:n+1),delta(1:l+1,1:m+1,1:n+1),epsilon(1:l+1,1:m+1,1:n+1)],...
                                [param_alpha(:,:,:,t+1), param_beta(:,:,:,t+1), param_gamma(:,:,:,t+1), param_delta(:,:,:,t+1), param_epsilon(:,:,:,t+1)]));

    disp('E = ')
    disp(E_value(t))

    if t > 1
        if E_value(t) > E_value(t-1)
            iteration = t;
            disp('iterationを終了します')
            break
        end
    end

end

l = 1;
m = 1;
n = 1;

f1_b(1) = 3 * param_alpha(1,1,1,iteration) + rho1 * (param_alpha(2,1,1,iteration) - param_alpha(1,1,1,iteration))...
                                           + rho2 * (param_alpha(1,2,1,iteration) - param_alpha(1,1,1,iteration))...
                                           + rho3 * (param_alpha(1,1,2,iteration) - param_alpha(1,1,1,iteration));

f2_b(1) = 3 * param_beta(1,1,1,iteration) + rho1 * (param_beta(2,1,1,iteration) - param_beta(1,1,1,iteration))...
                                          + rho2 * (param_beta(1,2,1,iteration) - param_beta(1,1,1,iteration))...
                                          + rho3 * (param_beta(1,1,2,iteration) - param_beta(1,1,1,iteration));  

f3_b(1) = 3 * param_gamma(1,1,1,iteration) + rho1 * (param_gamma(2,1,1,iteration) - param_gamma(1,1,1,iteration))...
                                           + rho2 * (param_gamma(1,2,1,iteration) - param_gamma(1,1,1,iteration))...
                                           + rho3 * (param_gamma(1,1,2,iteration) - param_gamma(1,1,1,iteration));

gmap_b(1) = 3 * param_delta(1,1,1,iteration) + rho1 * (param_delta(2,1,1,iteration) - param_delta(1,1,1,iteration))...
                                             + rho2 * (param_delta(1,2,1,iteration) - param_delta(1,1,1,iteration))...
                                             + rho3 * (param_delta(1,1,2,iteration) - param_delta(1,1,1,iteration));  

hmap_b(1) = 3 * param_epsilon(1,1,1,iteration) + rho1 * (param_epsilon(2,1,1,iteration) - param_epsilon(1,1,1,iteration))...
                                               + rho2 * (param_epsilon(1,2,1,iteration) - param_epsilon(1,1,1,iteration))...
                                               + rho3 * (param_epsilon(1,1,2,iteration) - param_epsilon(1,1,1,iteration));  


for j = 1:length(k) - 1

    l_now(j+1) = floor(si_b(j+1,1) / sigma1);
    rho1 = si_b(j+1,1) / sigma1  - l_now(j+1);

    m_now(j+1) = floor(si_b(j+1,2) / sigma2);
    rho2 = si_b(j+1,2) / sigma2  - m_now(j+1);

    n_now(j+1) = floor(si_b(j+1,3) / sigma3);
    rho3 = si_b(j+1,3) / sigma3  - n_now(j+1);

    if l_now(j+1) > l_now(j)
        l_num(j+1) = l_num(j) + 1; % 対応する格子点に番号をつけていく
    else
        l_num(j+1) = l_num(j);
    end

    if m_now(j+1) > m_now(j)
        m_num(j+1) = m_num(j) + 1; % 対応する格子点に番号をつけていく
    else
        m_num(j+1) = m_num(j);
    end

    if n_now(j+1) > n_now(j)
        n_num(j+1) = n_num(j) + 1; % 対応する格子点に番号をつけていく
    else
        n_num(j+1) = n_num(j);
    end

    l = l_num(j); %代入しやすくするため
    m = m_num(j);
    n = n_num(j);

    f1_b(j+1) = 3 * param_alpha(l,m,n,iteration) + rho1 * (param_alpha(l+1,m,n,iteration) - param_alpha(l,m,n,iteration))...
                                                 + rho2 * (param_alpha(l,m+1,n,iteration) - param_alpha(l,m,n,iteration))...
                                                 + rho3 * (param_alpha(l,m,n+1,iteration) - param_alpha(l,m,n,iteration));

    f2_b(j+1) = 3 * param_beta(l,m,n,iteration) + rho1 * (param_beta(l+1,m,n,iteration) - param_beta(l,m,n,iteration))...
                                                + rho2 * (param_beta(l,m+1,n,iteration) - param_beta(l,m,n,iteration))...
                                                + rho3 * (param_beta(l,m,n+1,iteration) - param_beta(l,m,n,iteration));  

    f3_b(j+1) = 3 * param_gamma(l,m,n,iteration) + rho1 * (param_gamma(l+1,m,n,iteration) - param_gamma(l,m,n,iteration))...
                                                 + rho2 * (param_gamma(l,m+1,n,iteration) - param_gamma(l,m,n,iteration))...
                                                 + rho3 * (param_gamma(l,m,n+1,iteration) - param_gamma(l,m,n,iteration));

    gmap_b(j+1) = 3 * param_delta(l,m,n,iteration) + rho1 * (param_delta(l+1,m,n,iteration) - param_delta(l,m,n,iteration))...
                                                   + rho2 * (param_delta(l,m+1,n,iteration) - param_delta(l,m,n,iteration))...
                                                   + rho3 * (param_delta(l,m,n+1,iteration) - param_delta(l,m,n,iteration));  

    hmap_b(j+1) = 3 * param_epsilon(l,m,n,iteration) + rho1 * (param_epsilon(l+1,m,n,iteration) - param_epsilon(l,m,n,iteration))...
                                                     + rho2 * (param_epsilon(l,m+1,n,iteration) - param_epsilon(l,m,n,iteration))...
                                                     + rho3 * (param_epsilon(l,m,n+1,iteration) - param_epsilon(l,m,n,iteration)); 


end



% 推定結果のplot--------------------------------------

figure;
hold on;
grid on;

axis([-1.7 1.7 -10 10]) % π/2 ≒ 1.57

plot(si_b(:,1), si_b(:,1), '--', si_b(:,1), f1_b(:),'LineWidth', 1.5) %z2 = f(s3) = tan(s3) の答え合わせ
xlabel('s1 = x')
ylabel('z1 = f1(s)')
legend('真値：s1','推定値：z1 = f1(s)')

hold off;


figure;
hold on;
grid on;

axis([-1.7 1.7 -10 10]) % π/2 ≒ 1.57

plot(si_b(:,3), si_b(:,3), '--', si_b(:,3), f2_b(:),'LineWidth', 1.5) %z2 = f(s3) = tan(s3) の答え合わせ
xlabel('s3 = θ')
ylabel('z2 = f2(s)')
legend('真値：tan(s3)','推定値：z2 = f2(s)')

hold off;


figure;
hold on;
grid on;

axis([-1.7 1.7 -10 10]) % π/2 ≒ 1.57

plot(si_b(:,2), si_b(:,2), '--', si_b(:,2), f3_b(:),'LineWidth', 1.5) %z2 = f(s3) = tan(s3) の答え合わせ
xlabel('s2 = y')
ylabel('z3 = f3(s)')
legend('真値：s2','推定値：z3 = f3(s)')

hold off;


figure;
hold on;
grid on;
axis([-1.7 1.7 0.5 10])

cos_3 = zeros(1,length(k));

for i=1:length(k)
    g_ans(i) = 1 / (cos(si_b(i,3)) * cos(si_b(i,3)) * cos(si_b(i,3)));
end

plot(si_b(:,3), g_ans(:), '--', si_b(:,3), gmap_b(:),'LineWidth', 1.5)
xlabel('s3 = θ')
ylabel('g(s)')
legend('真値：1/cos^3(s3)','推定値：g(s)')

hold off;


figure;
hold on;
grid on;
axis([-1.7 1.7 -0.2 1.2])

plot(si_b(:,3), cos(si_b(:,3)), '--', si_b(:,3), hmap_b(:),'LineWidth', 1.5)
xlabel('s3 = θ')
ylabel('h(s)')
legend('真値：cos(s3)','推定値：h(s)')

hold off;
