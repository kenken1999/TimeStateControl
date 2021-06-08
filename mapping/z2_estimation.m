clear;
close all;

dk = 0.05; %%時間刻み
Kfin = 1.2; %シミュレーション終了時間
k = [0:dk:Kfin];

u1 = ones(1,length(k)) * 2;
u2 = ones(1,length(k)) * 2;

si = zeros(length(k),3); %観測するセンサ変数 , 答えは(s1, s2, s3)=(x ,y, θ)
si(1,:) = [0 0 0]; %(s1, s2, s3)=(x ,y, θ)の初期値を設定

zi = zeros(length(k),3); %変換後の状態変数 (z1, z2, z3)=(x, tanθ, y)

p_now = zeros(1,length(k));

alpha = sym('alpha',[1 300]);

sigma = 0.05; %スケーリング定数
p = 1;
E = 0;

for j = 1:length(k)-1

    si(j+1,3) = si(j,3) + u2(j) * dk;
    si(j+1,1) = si(j,1) + u1(j) * cos(si(j+1,3)) * dk;
    si(j+1,2) = si(j,2) + u1(j) * sin(si(j+1,3)) * dk;

    zi(j+1,1) = si(j+1,1); %z1=s1は既知
    zi(j+1,3) = si(j+1,2); %z3=s2は既知

    p_now(j+1) = floor(si(j+1,3) / sigma) % p = 時刻kのi
    u = si(j+1,3) / sigma  - p_now(j+1);

    if p_now(j+1) > p_now(j)
        p = p + 1;
    end

    %zi(j,2) = alpha(p) + u * (alpha(p+1) - alpha(p)); %z2=f(s3)

    E = E + (alpha(p) + u * (alpha(p+1) - alpha(p)) - (zi(j+1,3)- zi(j,3)) / (zi(j+1,1)- zi(j,1))) ^ 2;   

end

length(k)
p

%正則化項の追加
% for i = 2:p

%     E = E + (alpha(i+1) - 2 * alpha(i) + alpha(i-1)) ^ 2;

% end

eta = 0.05; %学習率
iteration = 10; %繰り返し回数（最大）

param = zeros(iteration,p+1);

E_value = zeros(1,iteration);

syms 'alpha%d' [1 p+1]


for t = 1:iteration

    for m = 1:p+1

        DE = diff(E,alpha(m));
        
        DE2 = subs(DE, alpha(1:p+1), param(t,:));     

        param(t+1,m) = param(t,m) - eta * double(DE2);

    end

    E_value(t) = double(subs(E, alpha(1:p+1), param(t+1,:)));

    disp('E = ')
    disp(E_value(t))

    if t > 1 && E_value(t) > E_value(t-1)
        iteration = t;
        disp('iterationを終了します')
        break
    end

end

p = 1;

for j = 1:length(k) - 1

    %p_now(j+1) = floor(si(j+1,3) / sigma); % p = 時刻kのi
    %u = si(j+1,3) / sigma  - p_now(j+1);

    if p_now(j+1) > p_now(j)
        p = p + 1;
    end

    zi(j+1,2) = param(iteration,p) + u * (param(iteration,p+1) - param(iteration,p)); %z2=f(s3)

end

hold on;
grid on;

axis([-0.1 2.5 -10 10])

plot(si(:,3), tan(si(:,3)), '--', si(:,3), zi(:,2),'LineWidth', 1.5) %z2 = f(s3) = tan(s3) の答え合わせ
xlabel('s3 = θ')
ylabel('z2 = f(s3)')
legend('真値：tan(s3)','推定値：z2 = f(s3)')

