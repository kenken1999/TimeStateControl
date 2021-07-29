clear;
close all;

tic

%all(f1,f2,f3,g,h)_estimation----------------------------------------------------

dk = 0.02;   %時間刻み
Kfin = 0.08; %シミュレーション終了時間
k = [0:dk:Kfin];

u1_b = ones(length(k),1) * 5;
u2_b = ones(length(k),1) * 5;

si_b = zeros(length(k),3); %観測するセンサ変数 , s = (s1, s2, s3) = (x ,y, θ)
si_b(1,:) = [0 0 0];   %(s1, s2, s3)=(x ,y, θ)の初期値を設定

f1_b = zeros(length(k),1);   %写像f1:s→z1の推定
f2_b = zeros(length(k),1);   %写像f2:s→z2の推定
f3_b = zeros(length(k),1);   %写像f3:s→z3の推定
gmap_b = zeros(length(k),1); %dz2/dz1 = μ2 = u2/u1 * g(s)の推定, g(s) = 1/cos(s3)^3
hmap_b = zeros(length(k),1); %dz1/dt = μ1 = u1*h(s)の推定, h(s) = cos(s3)

l_now = zeros(length(k),1);  % l_now = 時刻kのl, 値が飛び飛び or 被る可能性あり
m_now = zeros(length(k),1);  % m_now = 時刻kのl, 値が飛び飛び or 被る可能性あり
n_now = zeros(length(k),1);  % n_now = 時刻kのl, 値が飛び飛び or 被る可能性あり

l_num = ones(length(k),1); %推定した格子点に前から順に番号をつけていく
m_num = ones(length(k),1);
n_num = ones(length(k),1);

rho1 = zeros(length(k),1);
rho2 = zeros(length(k),1);
rho3 = zeros(length(k),1);


alpha = sym('alpha',[50 50 50]); %l,m,nの順
beta = sym('beta',[50 50 50]);
gamma = sym('gamma',[50 50 50]);
delta = sym('delta',[50 50 50]);
epsilon = sym('epsilon',[50 50 50]);

sigma1 = 0.15; %s1のスケーリング定数
sigma2 = 0.15; %s2のスケーリング定数
sigma3 = 0.15; %s3のスケーリング定数

E = 0;
E1 = 0;
E2 = 0;
E3 = 0;
E4 = 0;
E4_1 = 0;
E4_2 = 0;
E4_3 = 0;
E4_4 = 0;
E4_5 = 0;

l_now(1) = floor(si_b(1,1) / sigma1);
rho1(1) = si_b(1,1) / sigma1  - l_now(1);

m_now(1) = floor(si_b(1,2) / sigma2);
rho2(1) = si_b(1,2) / sigma2  - m_now(1);

n_now(1) = floor(si_b(1,3) / sigma3);
rho3(1) = si_b(1,3) / sigma3  - n_now(1);


for j = 1:length(k) - 1

    si_b(j+1,3) = si_b(j,3) + u2_b(j) * dk;
    si_b(j+1,1) = si_b(j,1) + u1_b(j) * cos(si_b(j+1,3)) * dk;
    si_b(j+1,2) = si_b(j,2) + u1_b(j) * sin(si_b(j+1,3)) * dk;

    l_now(j+1) = floor(si_b(j+1,1) / sigma1);
    rho1(j+1) = si_b(j+1,1) / sigma1  - l_now(j+1);

    m_now(j+1) = floor(si_b(j+1,2) / sigma2);
    rho2(j+1) = si_b(j+1,2) / sigma2  - m_now(j+1);

    n_now(j+1) = floor(si_b(j+1,3) / sigma3);
    rho3(j+1) = si_b(j+1,3) / sigma3  - n_now(j+1);

    if l_now(j+1) == l_now(j)
        l_num(j+1) = l_num(j);
    else
        l_num(j+1) = l_num(j) + 1; % 対応する格子点に番号をつけていく
    end

    if m_now(j+1) == m_now(j) 
        m_num(j+1) = m_num(j);
    else
        m_num(j+1) = m_num(j) + 1; % 対応する格子点に番号をつけていく
    end

    if n_now(j+1) == n_now(j)
        n_num(j+1) = n_num(j);
    else
        n_num(j+1) = n_num(j) + 1; % 対応する格子点に番号をつけていく
    end

    l = l_num(j); %代入しやすくするため
    m = m_num(j);
    n = n_num(j);

    l2 = l_num(j+1); %代入しやすくするため
    m2 = m_num(j+1);
    n2 = n_num(j+1);

    E1 = E1 + ( beta(l,m,n) + rho1(j) * (beta(l+1,m,n) - beta(l,m,n)) + rho2(j) * (beta(l,m+1,n) - beta(l,m,n)) + rho3(j) * (beta(l,m,n+1) - beta(l,m,n))...
                 - ( (gamma(l2,m2,n2) + rho1(j+1) * (gamma(l2+1,m2,n2) - gamma(l2,m2,n2)) + rho2(j+1) * (gamma(l2,m2+1,n2) - gamma(l2,m2,n2)) + rho3(j+1) * (gamma(l2,m2,n2+1) - gamma(l2,m2,n2)))...
                 - (gamma(l,m,n) + rho1(j) * (gamma(l+1,m,n) - gamma(l,m,n)) + rho2(j) * (gamma(l,m+1,n) - gamma(l,m,n)) + rho3(j) * (gamma(l,m,n+1) - gamma(l,m,n))) )...
                 / ( (alpha(l2,m2,n2) + rho1(j+1) * (alpha(l2+1,m2,n2) - alpha(l2,m2,n2)) + rho2(j+1) * (alpha(l2,m2+1,n2) - alpha(l2,m2,n2)) + rho3(j+1) * (alpha(l2,m2,n2+1) - alpha(l2,m2,n2)))...
                 - (alpha(l,m,n) + rho1(j) * (alpha(l+1,m,n) - alpha(l,m,n)) + rho2(j) * (alpha(l,m+1,n) - alpha(l,m,n)) + rho3(j) * (alpha(l,m,n+1) - alpha(l,m,n))) ) ) ^ 2;

    E2 = E2 + ( u2_b(j) / u1_b(j) * ( delta(l,m,n) + rho1(j) * (delta(l+1,m,n) - delta(l,m,n)) + rho2(j) * (delta(l,m+1,n) - delta(l,m,n)) + rho3(j) * (delta(l,m,n+1) - delta(l,m,n)) )...
                 - ( (beta(l2,m2,n2) + rho1(j+1) * (beta(l2+1,m2,n2) - beta(l2,m2,n2)) + rho2(j+1) * (beta(l2,m2+1,n2) - beta(l2,m2,n2)) + rho3(j+1) * (beta(l2,m2,n2+1) - beta(l2,m2,n2)))...
                 - (beta(l,m,n) + rho1(j) * (beta(l+1,m,n) - beta(l,m,n)) + rho2(j) * (beta(l,m+1,n) - beta(l,m,n)) + rho3(j) * (beta(l,m,n+1) - beta(l,m,n))) )...
                 / ( (alpha(l2,m2,n2) + rho1(j+1) * (alpha(l2+1,m2,n2) - alpha(l2,m2,n2)) + rho2(j+1) * (alpha(l2,m2+1,n2) - alpha(l2,m2,n2)) + rho3(j+1) * (alpha(l2,m2,n2+1) - alpha(l2,m2,n2)))...
                 - (alpha(l,m,n) + rho1(j) * (alpha(l+1,m,n) - alpha(l,m,n)) + rho2(j) * (alpha(l,m+1,n) - alpha(l,m,n)) + rho3(j) * (alpha(l,m,n+1) - alpha(l,m,n))) ) ) ^ 2;
    
    E3 = E3 + ( u1_b(j) * ( epsilon(l,m,n) + rho1(j) * (epsilon(l+1,m,n) - epsilon(l,m,n)) + rho2(j) * (epsilon(l,m+1,n) - epsilon(l,m,n)) + rho3(j) * (epsilon(l,m,n+1) - epsilon(l,m,n)) )...
                 - ( (alpha(l2,m2,n2) + rho1(j+1) * (alpha(l2+1,m2,n2) - alpha(l2,m2,n2)) + rho2(j+1) * (alpha(l2,m2+1,n2) - alpha(l2,m2,n2)) + rho3(j+1) * (alpha(l2,m2,n2+1) - alpha(l2,m2,n2)))...
                 - (alpha(l,m,n) + rho1(j) * (alpha(l+1,m,n) - alpha(l,m,n)) + rho2(j) * (alpha(l,m+1,n) - alpha(l,m,n)) + rho3(j) * (alpha(l,m,n+1) - alpha(l,m,n))) )...
                 / dk ) ^ 2;

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

disp('l2 = ')
disp(l2)

disp('m2 = ')
disp(m2)

disp('n2 = ')
disp(n2)



%---正則化項の追加---------------

xi_1 = 1;
xi_2 = 1;
xi_3 = 1;
xi_4 = 1;
xi_5 = 1;


for a = 1 : l2 - 1
    for b = 1 : m2 - 1
        for c = 1 : n2 - 1

            E4_1 = (alpha(a+2,b,c) - 2 * alpha(a+1,b,c) - alpha(a,b,c)) ^ 2 + (alpha(a,b+2,c) - 2 * alpha(a,b+1,c) - alpha(a,b,c)) ^ 2 + (alpha(a,b,c+2) - 2 * alpha(a,b,c+1) - alpha(a,b,c)) ^ 2;
            E4_2 = (beta(a+2,b,c) - 2 * beta(a+1,b,c) + beta(a,b,c)) ^ 2 + (beta(a,b+2,c) - 2 * beta(a,b+1,c) + beta(a,b,c)) ^ 2 + (beta(a,b,c+2) - 2 * beta(a,b,c+1) + beta(a,b,c)) ^ 2;
            E4_3 = (gamma(a+2,b,c) - 2 * gamma(a+1,b,c) + gamma(a,b,c)) ^ 2 + (gamma(a,b+2,c) - 2 * gamma(a,b+1,c) + gamma(a,b,c)) ^ 2 + (gamma(a,b,c+2) - 2 * gamma(a,b,c+1) + gamma(a,b,c)) ^ 2;
            E4_4 = (delta(a+2,b,c) - 2 * delta(a+1,b,c) + delta(a,b,c)) ^ 2 + (delta(a,b+2,c) - 2 * delta(a,b+1,c) + delta(a,b,c)) ^ 2 + (delta(a,b,c+2) - 2 * delta(a,b,c+1) + delta(a,b,c)) ^ 2;
            E4_5 = (epsilon(a+2,b,c) - 2 * epsilon(a+1,b,c) + epsilon(a,b,c)) ^ 2 + (epsilon(a,b+2,c) - 2 * epsilon(a,b+1,c) + epsilon(a,b,c)) ^ 2 + (epsilon(a,b,c+2) - 2 * epsilon(a,b,c+1) + epsilon(a,b,c)) ^ 2;

            % E4_1 = (alpha(a+1,b,c) - 2 * alpha(a,b,c) + alpha(a-1,b,c)) ^ 2 + (alpha(a,b+1,c) - 2 * alpha(a,b,c) + alpha(a,b-1,c)) ^ 2 + (alpha(a,b,c+1) - 2 * alpha(a,b,c) + alpha(a,b,c-1)) ^ 2;
            % E4_2 = (beta(a+1,b,c) - 2 * beta(a,b,c) + beta(a-1,b,c)) ^ 2 + (beta(a,b+1,c) - 2 * beta(a,b,c) + beta(a,b-1,c)) ^ 2 + (beta(a,b,c+1) - 2 * beta(a,b,c) + beta(a,b,c-1)) ^ 2;
            % E4_3 = (gamma(a+1,b,c) - 2 * gamma(a,b,c) + gamma(a-1,b,c)) ^ 2 + (gamma(a,b+1,c) - 2 * gamma(a,b,c) + gamma(a,b-1,c)) ^ 2 + (gamma(a,b,c+1) - 2 * gamma(a,b,c) + gamma(a,b,c-1)) ^ 2;
            % E4_4 = (delta(a+1,b,c) - 2 * delta(a,b,c) + delta(a-1,b,c)) ^ 2 + (delta(a,b+1,c) - 2 * delta(a,b,c) + delta(a,b-1,c)) ^ 2 + (delta(a,b,c+1) - 2 * delta(a,b,c) + delta(a,b,c-1)) ^ 2;
            % E4_5 = (epsilon(a+1,b,c) - 2 * epsilon(a,b,c) + epsilon(a-1,b,c)) ^ 2 + (epsilon(a,b+1,c) - 2 * epsilon(a,b,c) + epsilon(a,b-1,c)) ^ 2 + (epsilon(a,b,c+1) - 2 * epsilon(a,b,c) + epsilon(a,b,c-1)) ^ 2;

            E4 = E4 + xi_1 * E4_1 + xi_2 * E4_2 + xi_3 * E4_3 + xi_4 * E4_4 + xi_5 * E4_5;


        end
    end
end


%---写像 f1,f2,f3,g,h の推定----------------------------

eta_f1 = 1.0 * 10 ^ (-5); %学習率 5
eta_f2 = 1.0 * 10 ^ (-4); %学習率 4
eta_f3 = 1.0 * 10 ^ (-5); %学習率 5
eta_g = 1.0 * 10 ^ (-5); %学習率 5
eta_h = 1.0 * 10 ^ (-4); %学習率 4

iteration = 50; %パラメータ更新回数（最大）

% param_alpha = rand([l2+1,m2+1,n2+1,iteration+1]);
% param_beta = rand([l2+1,m2+1,n2+1,iteration+1]);
% param_gamma = rand([l2+1,m2+1,n2+1,iteration+1]);
% param_delta = rand([l2+1,m2+1,n2+1,iteration+1]);
% param_epsilon = rand([l2+1,m2+1,n2+1,iteration+1]);

param_alpha = rand([l2+1,m2+1,n2+1,iteration+1]) * 0.1 - 0.05;
param_beta = rand([l2+1,m2+1,n2+1,iteration+1]) * 0.1 - 0.05;
param_gamma = rand([l2+1,m2+1,n2+1,iteration+1]) * 0.1 - 0.05;
param_delta = rand([l2+1,m2+1,n2+1,iteration+1]) * 0.1 - 0.05;
param_epsilon = rand([l2+1,m2+1,n2+1,iteration+1]) * 0.1 - 0.05;

%---Eの設定-----------------------------

E1_initial = double(subs(E1, [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                             [param_alpha(:,:,:,1), param_beta(:,:,:,1), param_gamma(:,:,:,1), param_delta(:,:,:,1), param_epsilon(:,:,:,1)]));

E2_initial = double(subs(E2, [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                             [param_alpha(:,:,:,1), param_beta(:,:,:,1), param_gamma(:,:,:,1), param_delta(:,:,:,1), param_epsilon(:,:,:,1)]));

E3_initial = double(subs(E3, [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                             [param_alpha(:,:,:,1), param_beta(:,:,:,1), param_gamma(:,:,:,1), param_delta(:,:,:,1), param_epsilon(:,:,:,1)]));

E4_initial = double(subs(E4, [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                             [param_alpha(:,:,:,1), param_beta(:,:,:,1), param_gamma(:,:,:,1), param_delta(:,:,:,1), param_epsilon(:,:,:,1)]));
                             

% zeta1 = 1; %E1の調整係数
% zeta2 = 0.1; %E2の調整係数
% zeta3 = 1; %E3の調整係数

zeta1 = 1 / (E1_initial + 1); %E1の調整係数
zeta2 = 1 / (E2_initial + 1); %E2の調整係数
zeta3 = 1 / (E3_initial + 1); %E3の調整係数

xi_all =  1 / (E4_initial + 1);

E = zeta1 * E1 + zeta2 * E2 + zeta3 * E3 + xi_all * E4;

E_initial = double(subs(E, [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                             [param_alpha(:,:,:,1), param_beta(:,:,:,1), param_gamma(:,:,:,1), param_delta(:,:,:,1), param_epsilon(:,:,:,1)]));


disp('E1 = ')
disp(E1_initial)
disp('zeta1 = ')
disp(zeta1)
disp('E2 = ')
disp(E2_initial)
disp('zeta2 = ')
disp(zeta2)
disp('E3 = ')
disp(E3_initial)
disp('zeta3 = ')
disp(zeta3)
disp('E4 = ')
disp(E4_initial)
disp('xi_all = ')
disp(xi_all)
disp('E = ')
disp(E_initial)
disp('--------------------')



E_value = zeros(1,iteration);
E1_value = zeros(1,iteration);
E2_value = zeros(1,iteration);
E3_value = zeros(1,iteration);
E4_value = zeros(1,iteration);

DEf1 = sym('DEf1',[l+1 m+1 n+1]);
DEf2 = sym('DEf2',[l+1 m+1 n+1]);
DEf3 = sym('DEf3',[l+1 m+1 n+1]);
DEg = sym('DEg',[l+1 m+1 n+1]);
DEh = sym('DEh',[l+1 m+1 n+1]);


for a = 1:l2+1
    for b = 1:m2+1
        for c = 1:n2+1

                DEf1(a,b,c) = diff(E,alpha(a,b,c));
                DEf2(a,b,c) = diff(E,beta(a,b,c));
                DEf3(a,b,c) = diff(E,gamma(a,b,c));
                DEg(a,b,c) = diff(E,delta(a,b,c));
                DEh(a,b,c) = diff(E,epsilon(a,b,c));

                DEf1_1{a,b,c} = matlabFunction(DEf1(a,b,c), 'vars', {alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)});
                DEf2_1{a,b,c} = matlabFunction(DEf2(a,b,c), 'vars', {alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)});
                DEf3_1{a,b,c} = matlabFunction(DEf3(a,b,c), 'vars', {alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)});
                DEg_1{a,b,c} = matlabFunction(DEg(a,b,c), 'vars', {alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)});
                DEh_1{a,b,c} = matlabFunction(DEh(a,b,c), 'vars', {alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)});

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

                DEf1_2 = DEf1_1{a,b,c}(param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t));
                DEf2_2 = DEf2_1{a,b,c}(param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t));
                DEf3_2 = DEf3_1{a,b,c}(param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t));
                DEg_2 = DEg_1{a,b,c}(param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t));
                DEh_2 = DEh_1{a,b,c}(param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t));

                % DEf1_2 = subs(DEf1(a,b,c), [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                %                            [param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t)]);

                % DEf2_2 = subs(DEf2(a,b,c), [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                %                            [param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t)]);

                % DEf3_2 = subs(DEf3(a,b,c), [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                %                            [param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t)]);

                % DEg_2 = subs(DEg(a,b,c), [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                %                          [param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t)]);

                % DEh_2 = subs(DEh(a,b,c), [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                %                          [param_alpha(:,:,:,t), param_beta(:,:,:,t), param_gamma(:,:,:,t), param_delta(:,:,:,t), param_epsilon(:,:,:,t)]);


                param_alpha(a,b,c,t+1) = param_alpha(a,b,c,t) - eta_f1 * DEf1_2;
                param_beta(a,b,c,t+1) = param_beta(a,b,c,t) - eta_f2 * DEf2_2;
                param_gamma(a,b,c,t+1) = param_gamma(a,b,c,t) - eta_f3 * DEf3_2;
                param_delta(a,b,c,t+1) = param_delta(a,b,c,t) - eta_g * DEg_2;
                param_epsilon(a,b,c,t+1) = param_epsilon(a,b,c,t) - eta_h * DEh_2;

               
            end
        end
    end

    param_alpha(1,1,1,t+1) = 0;
    param_beta(1,1,1,t+1) = 0;
    param_gamma(1,1,1,t+1) = 0;
    param_delta(1,1,1,t+1) = 1;
    param_epsilon(1,1,1,t+1) = 1;


    E1_value(t) = double(subs(E1, [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                                 [param_alpha(:,:,:,t+1), param_beta(:,:,:,t+1), param_gamma(:,:,:,t+1), param_delta(:,:,:,t+1), param_epsilon(:,:,:,t+1)]));

    E2_value(t) = double(subs(E2, [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                                 [param_alpha(:,:,:,t+1), param_beta(:,:,:,t+1), param_gamma(:,:,:,t+1), param_delta(:,:,:,t+1), param_epsilon(:,:,:,t+1)]));

    E3_value(t) = double(subs(E3, [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                                 [param_alpha(:,:,:,t+1), param_beta(:,:,:,t+1), param_gamma(:,:,:,t+1), param_delta(:,:,:,t+1), param_epsilon(:,:,:,t+1)]));
    
    E4_value(t) = double(subs(E4, [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                                 [param_alpha(:,:,:,t+1), param_beta(:,:,:,t+1), param_gamma(:,:,:,t+1), param_delta(:,:,:,t+1), param_epsilon(:,:,:,t+1)]));

    E_value(t) = double(subs(E, [alpha(1:l2+1,1:m2+1,1:n2+1),beta(1:l2+1,1:m2+1,1:n2+1),gamma(1:l2+1,1:m2+1,1:n2+1),delta(1:l2+1,1:m2+1,1:n2+1),epsilon(1:l2+1,1:m2+1,1:n2+1)],...
                                [param_alpha(:,:,:,t+1), param_beta(:,:,:,t+1), param_gamma(:,:,:,t+1), param_delta(:,:,:,t+1), param_epsilon(:,:,:,t+1)]));


    disp('t = ')
    disp(t)
    disp('E1 = ')
    disp(E1_value(t))
    disp('E2 = ')
    disp(E2_value(t))
    disp('E3 = ')
    disp(E3_value(t))
    disp('E4 = ')
    disp(E4_value(t))
    disp('E = ')
    disp(E_value(t))

    if t > 1
        if (E1_value(t) > E1_value(t-1))
            disp('E1が増加しました')
        end
        if (E2_value(t) > E2_value(t-1))
            disp('E2が増加しました')
        end
        if (E3_value(t) > E3_value(t-1))
            disp('E3が増加しました')
        end
        if (E4_value(t) > E4_value(t-1))
            disp('E4が増加しました')
        end
        if (E_value(t) > E_value(t-1))
            disp('Eが増加しました')
        end
        if (E1_value(t) > E1_value(t-1)) && (E2_value(t) > E2_value(t-1)) && (E3_value(t) > E3_value(t-1)) && (E4_value(t) > E4_value(t-1))
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



% f1_b(1) = 3 * param_alpha(1,1,1,iteration) - (1-rho1(1)) * (param_alpha(2,1,1,iteration) - param_alpha(1,1,1,iteration))...
%                                            - (1-rho2(1)) * (param_alpha(1,2,1,iteration) - param_alpha(1,1,1,iteration))...
%                                            - (1-rho3(1)) * (param_alpha(1,1,2,iteration) - param_alpha(1,1,1,iteration));

% f2_b(1) = 3 * param_beta(1,1,1,iteration) - (1-rho1(1)) * (param_beta(2,1,1,iteration) - param_beta(1,1,1,iteration))...
%                                           - (1-rho2(1)) * (param_beta(1,2,1,iteration) - param_beta(1,1,1,iteration))...
%                                           - (1-rho3(1)) * (param_beta(1,1,2,iteration) - param_beta(1,1,1,iteration));  

% f3_b(1) = 3 * param_gamma(1,1,1,iteration) - (1-rho1(1)) * (param_gamma(2,1,1,iteration) - param_gamma(1,1,1,iteration))...
%                                            - (1-rho2(1)) * (param_gamma(1,2,1,iteration) - param_gamma(1,1,1,iteration))...
%                                            - (1-rho3(1)) * (param_gamma(1,1,2,iteration) - param_gamma(1,1,1,iteration));

% gmap_b(1) = 3 * param_delta(1,1,1,iteration) - (1-rho1(1)) * (param_delta(2,1,1,iteration) - param_delta(1,1,1,iteration))...
%                                              - (1-rho2(1)) * (param_delta(1,2,1,iteration) - param_delta(1,1,1,iteration))...
%                                              - (1-rho3(1)) * (param_delta(1,1,2,iteration) - param_delta(1,1,1,iteration));  

% hmap_b(1) = 3 * param_epsilon(1,1,1,iteration) - (1-rho1(1)) * (param_epsilon(2,1,1,iteration) - param_epsilon(1,1,1,iteration))...
%                                                - (1-rho2(1)) * (param_epsilon(1,2,1,iteration) - param_epsilon(1,1,1,iteration))...
%                                                - (1-rho3(1)) * (param_epsilon(1,1,2,iteration) - param_epsilon(1,1,1,iteration));  


for j = 1:length(k)

    l = l_num(j); %代入しやすくするため
    m = m_num(j);
    n = n_num(j);

    f1_b(j) = param_alpha(l,m,n,iteration) + rho1(j) * (param_alpha(l+1,m,n,iteration) - param_alpha(l,m,n,iteration))...
                                                 + rho2(j) * (param_alpha(l,m+1,n,iteration) - param_alpha(l,m,n,iteration))...
                                                 + rho3(j) * (param_alpha(l,m,n+1,iteration) - param_alpha(l,m,n,iteration));

    f2_b(j) = param_beta(l,m,n,iteration) + rho1(j) * (param_beta(l+1,m,n,iteration) - param_beta(l,m,n,iteration))...
                                                + rho2(j) * (param_beta(l,m+1,n,iteration) - param_beta(l,m,n,iteration))...
                                                + rho3(j) * (param_beta(l,m,n+1,iteration) - param_beta(l,m,n,iteration));  

    f3_b(j) = param_gamma(l,m,n,iteration) + rho1(j) * (param_gamma(l+1,m,n,iteration) - param_gamma(l,m,n,iteration))...
                                                 + rho2(j) * (param_gamma(l,m+1,n,iteration) - param_gamma(l,m,n,iteration))...
                                                 + rho3(j) * (param_gamma(l,m,n+1,iteration) - param_gamma(l,m,n,iteration));

    gmap_b(j) = param_delta(l,m,n,iteration) + rho1(j) * (param_delta(l+1,m,n,iteration) - param_delta(l,m,n,iteration))...
                                                   + rho2(j) * (param_delta(l,m+1,n,iteration) - param_delta(l,m,n,iteration))...
                                                   + rho3(j) * (param_delta(l,m,n+1,iteration) - param_delta(l,m,n,iteration));  

    hmap_b(j) = param_epsilon(l,m,n,iteration) + rho1(j) * (param_epsilon(l+1,m,n,iteration) - param_epsilon(l,m,n,iteration))...
                                                     + rho2(j) * (param_epsilon(l,m+1,n,iteration) - param_epsilon(l,m,n,iteration))...
                                                     + rho3(j) * (param_epsilon(l,m,n+1,iteration) - param_epsilon(l,m,n,iteration)); 


end


toc

% 推定結果のplot--------------------------------------

figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_b(:,1), si_b(:,1), '--m', si_b(:,1), f1_b(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z1 = f1(s) = s1 の答え合わせ
xlabel('s1 = x')
ylabel('z1 = f1(s)')
legend('真値：s1','推定値：z1 = f1(s)')

hold off;


figure;
hold on;
grid on;

axis([-1.7 1.7 -10 10]) % π/2 ≒ 1.57

plot(si_b(:,3), tan(si_b(:,3)), '--m', si_b(:,3), f2_b(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z2 = f2(s) = tan(s3) の答え合わせ
xlabel('s3 = θ')
ylabel('z2 = f2(s)')
legend('真値：tan(s3)','推定値：z2 = f2(s)')

hold off;


figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_b(:,2), si_b(:,2), '--m', si_b(:,2), f3_b(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z3 = f3(s) = s2 の答え合わせ
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

plot(si_b(:,3), g_ans(:), '--m', si_b(:,3), gmap_b(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel('s3 = θ')
ylabel('g(s)')
legend('真値：1/cos^3(s3)','推定値：g(s)')

hold off;


figure;
hold on;
grid on;
axis([-1.7 1.7 -0.2 1.2])

plot(si_b(:,3), cos(si_b(:,3)), '--m', si_b(:,3), hmap_b(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5)
xlabel('s3 = θ')
ylabel('h(s)')
legend('真値：cos(s3)','推定値：h(s)')

hold off;
