clear;
close all;

dk = 0.1; %%時間刻み
Kfin = 5; %シミュレーション終了時間
k = [0:dk:Kfin];

u1 = ones(1,length(k)) * 5;
u2 = ones(1,length(k)) * 5;

si = zeros(length(k),3); %観測するセンサ変数 , 答えは(s1, s2, s3)=(x ,y, θ)
si(1,:) = [0 0 0]; %(s1, s2, s3)=(x ,y, θ)の初期値を設定

zi = zeros(length(k),3); %変換後の状態変数 (z1, z2, z3)=(x, tanθ, y)

p_now = zeros(1,length(k));

alpha = sym('alpha',[1 length(k)])

sigma = 0.01; %スケーリング定数
p = 0;
E = 0;

for j = 1:length(k)-1

    si(j+1,3) = si(j,3) + u2(j) * dk;
    si(j+1,1) = si(j,1) + u1(j) * cos(si(j+1,3)) * dk;
    si(j+1,2) = si(j,2) + u1(j) * sin(si(j+1,3)) * dk;

    zi(j+1,1) = si(j+1,1); %z1=s1は既知
    zi(j+1,3) = si(j+1,2); %z3=s2は既知

    p_now(j+1) = floor(si(j+1,3) / sigma); % p = 時刻kのi
    u = si(j+1,3) / sigma  - p_now(j+1);

    if p_now(j+1) > p_now(j)
        p = p + 1;
    end

    %zi(j,2) = alpha(p) + u * (alpha(p+1) - alpha(p)); %z2=f(s3)

    E = E + (alpha(p) + u * (alpha(p+1) - alpha(p)) - (zi(j+1,3)- zi(j,3)) / (zi(j+1,1)- zi(j,1))) ^ 2;

    

end

E

for i = 2:p

    E = E + (alpha(i+1) - 2 * alpha(i) + alpha(i-1)) ^ 2;

end

E

eta = 0.01; %学習率
iteration = 3; %繰り返し回数（最大）

param = ones(iteration,p);

for t=1:iteration-1

    for m = 1:p-1

        syms E

        Df = diff(E,alpha(m))

        Df2 = Df(param(t,m));
        param(t+1,m) = param(t,m) - eta * double(Df2);

    end

end

param(iteration,:)


