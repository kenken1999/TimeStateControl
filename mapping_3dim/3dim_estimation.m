clear;
close all;

%all(f1,f2,f3,g,h)_estimation----------------------------------------------------

dk = 0.02;   %時間刻み
Kfin = 0.62; %シミュレーション終了時間
k = [0:dk:Kfin];

u1_b = ones(1,length(k)) * 5;
u2_b = ones(1,length(k)) * 5;

si_b = zeros(length(k),3); %観測するセンサ変数 , s = (s1, s2, s3) = (x ,y, θ)
si_b(1,:) = [0 0 -pi/2];   %(s1, s2, s3)=(x ,y, θ)の初期値を設定

f1_b = zeros(length(k),3);   %写像f1:s→z1の推定
f2_b = zeros(length(k),3);   %写像f2:s→z2の推定
f3_b = zeros(length(k),3);   %写像f3:s→z3の推定
g_b = zeros(1,length(k)); %dz2/dz1 = μ2 = u2/u1 * g(s)の推定, g(s) = 1/cos(s3)^3
h_b = zeros(1,length(k)); %dz1/dt = μ1 = u1*h(s)の推定, h(s) = cos(s3)

l_now = zeros(1,length(k));  % l_now = 時刻kのl, 値が飛び飛び or 被る可能性あり
m_now = zeros(1,length(k));  % m_now = 時刻kのl, 値が飛び飛び or 被る可能性あり
n_now = zeros(1,length(k));  % n_now = 時刻kのl, 値が飛び飛び or 被る可能性あり

l = 1; %推定した格子点に前から順に番号をつけていく
m = 1;
n = 1;

alpha = sym('alpha',[1 300]);
beta = sym('beta',[1 300]);
gamma = sym('gamma',[1 300]);
delta = sym('delta',[1 300]);
epsilon = sym('epsilon',[1 300]);

sigma1 = 0.1; %s1のスケーリング定数
sigma2 = 0.1; %s2のスケーリング定数
sigma3 = 0.1; %s3のスケーリング定数

E = 0;

for j = 1:length(k) - 1

    si_b(j+1,3) = si_b(j,3) + u2_b(j) * dk;
    si_b(j+1,1) = si_b(j,1) + u1_b(j) * cos(si_b(j+1,3)) * dk;
    si_b(j+1,2) = si_b(j,2) + u1_b(j) * sin(si_b(j+1,3)) * dk;

    zi_b(j+1,1) = si_b(j+1,1); %z1=s1は既知
    zi_b(j+1,3) = si_b(j+1,2); %z3=s2は既知

    p_now(j+1) = floor(si_b(j+1,3) / sigma); % p_now = 時刻kのi, 値が飛び飛び or 被る可能性あり
    u = si_b(j+1,3) / sigma  - p_now(j+1);

    if p_now(j+1) > p_now(j)
        p = p + 1; % 対応する格子点に番号をつけていく
    end

    Ef = Ef + ( alpha(p) + u * (alpha(p+1) - alpha(p)) - (zi_b(j+1,3) - zi_b(j,3)) / (zi_b(j+1,1) - zi_b(j,1)) ) ^ 2;

    Eh = Eh + ( u1_b(j) * (gamma(p) + u * (gamma(p+1) - gamma(p))) - (zi_b(j+1,1) - zi_b(j,1)) / dk )^ 2;

end

p_now


%---正則化項の追加---------------------------------
% for i = 2:p

%     Ef = Ef + (alpha(i+1) - 2 * alpha(i) + alpha(i-1)) ^ 2;

% end


%---写像 f, g, hの推定----------------------------

eta_f = 0.05; %学習率
eta_g = 0.5; %学習率
eta_h = 0.01; %学習率

iteration = 30; %パラメータ更新回数（最大）

param_alpha = zeros(iteration,p+1);
param_beta = zeros(iteration,p+1);
param_gamma = zeros(iteration,p+1);

Ef_value = zeros(1,iteration);
Eh_value = zeros(1,iteration);

syms 'alpha%d' [1 p+1]
syms 'beta%d' [1 p+1]
syms 'gamma%d' [1 p+1]


%---写像 f, hの推定-----------------

for t = 1:iteration

    for m = 1:p+1

        DEf = diff(Ef,alpha(m));
        DEh = diff(Eh,gamma(m));
        
        DEf2 = subs(DEf, alpha(1:p+1), param_alpha(t,:));
        DEh2 = subs(DEh, gamma(1:p+1), param_gamma(t,:));

        param_alpha(t+1,m) = param_alpha(t,m) - eta_f * double(DEf2);
        param_gamma(t+1,m) = param_gamma(t,m) - eta_h * double(DEh2);

    end

    Ef_value(t) = double(subs(Ef, alpha(1:p+1), param_alpha(t+1,:)));
    Eh_value(t) = double(subs(Eh, gamma(1:p+1), param_gamma(t+1,:)));

    disp('Ef = ')
    disp(Ef_value(t))

    disp('Eh = ')
    disp(Eh_value(t))

    if t > 1
        if Ef_value(t) > Ef_value(t-1) || Eh_value(t) > Eh_value(t-1)
            iteration = t;
            disp('iterationを終了します')
            break
        end
    end

end

p = 1;

zi_b(1,2) = param_alpha(iteration,1) + u * (param_alpha(iteration,2) - param_alpha(iteration,1)); %z2=f(s3)
hmap_b(1) = param_gamma(iteration,1) + u * (param_gamma(iteration,2) - param_gamma(iteration,1)); %h(s3)

for j = 1:length(k) - 1

    p_now(j+1) = floor(si_b(j+1,3) / sigma);% p_now = 時刻kのi, 値が飛び飛び or 被る可能性あり
    u = si_b(j+1,3) / sigma  - p_now(j+1);

    if p_now(j+1) > p_now(j)
        p = p + 1; 
    end

    zi_b(j+1,2) = param_alpha(iteration,p) + u * (param_alpha(iteration,p+1) - param_alpha(iteration,p)); %z2=f(s3)
    hmap_b(j+1) = param_gamma(iteration,p) + u * (param_gamma(iteration,p+1) - param_gamma(iteration,p)); %h(s3)

    %推定したz2_bを利用して， Egを生成
    Eg = Eg + ( u2_b(j) / u1_b(j) * (beta(p) + u * (beta(p+1) - beta(p))) - (zi_b(j+1,2)- zi_b(j,2)) / (zi_b(j+1,1)- zi_b(j,1)) ) ^ 2;

end


%---写像gの推定-----------------

for t = 1:iteration

    for m = 1:p+1

        DEg = diff(Eg,beta(m));
        
        DEg2 = subs(DEg, beta(1:p+1), param_beta(t,:));

        param_beta(t+1,m) = param_beta(t,m) - eta_g * double(DEg2);

    end

    Eg_value(t) = double(subs(Eg, beta(1:p+1), param_beta(t+1,:)));

    disp('Eg = ')
    disp(Eg_value(t))

    if t > 1
        if Eg_value(t) > Eg_value(t-1)
            iteration = t;
            disp('iterationを終了します')
            break
        end
    end

end


p = 1;

gmap_b(1) = param_beta(iteration,1) + u * (param_beta(iteration,2) - param_beta(iteration,1)); %h(s3)


for j = 1:length(k) - 1

    p_now(j+1) = floor(si_b(j+1,3) / sigma);% p_now = 時刻kのi, 値が飛び飛び or 被る可能性あり
    u = si_b(j+1,3) / sigma  - p_now(j+1);

    if p_now(j+1) > p_now(j)
        p = p + 1; 
    end

    gmap_b(j+1) = param_beta(iteration,p) + u * (param_beta(iteration,p+1) - param_beta(iteration,p)); %g(s3)

end


% 推定結果のplot--------------------------------------

figure;
hold on;
grid on;

axis([-1.7 1.7 -10 10]) % π/2 ≒ 1.57

plot(si_b(:,3), tan(si_b(:,3)), '--', si_b(:,3), zi_b(:,2),'LineWidth', 1.5) %z2 = f(s3) = tan(s3) の答え合わせ
xlabel('s3 = θ')
ylabel('z2 = f(s3)')
legend('真値：tan(s3)','推定値：z2 = f(s3)')

hold off;


figure;
hold on;
grid on;
axis([-1.7 1.7 -0.2 1.2])

plot(si_b(:,3), cos(si_b(:,3)), '--', si_b(:,3), hmap_b(:),'LineWidth', 1.5)
xlabel('s3 = θ')
ylabel('h(s3)')
legend('真値：cos(s3)','推定値：h(s3)')

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
ylabel('g(s3)')
legend('真値：1/cos^3(s3)','推定値：g(s3)')

hold off;

