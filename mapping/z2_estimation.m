clear;
close all;

dk = 0.1; %%時間刻み
Kfin = 10; %シミュレーション終了時間
k = [0:dk:Kfin];

u1 = ones(1,length(k)) * 5;
u2 = ones(1,length(k)) * 5;

si = zeros(length(k),3); %観測するセンサ変数 , 答えは(s1, s2, s3)=(x ,y, θ)
si(1,:) = [0 0 0]; %(s1, s2, s3)=(x ,y, θ)の初期値を設定

zi = zeros(length(k),3); %変換後の状態変数 (z1, z2, z3)=(x, tanθ, y)


for j = 1:length(k)-1

    zi(j,1) = si(j,1); %z1=s1は既知
    zi(j,3) = si(j,2); %z3=s2は既知

    sigma = 0.01; %スケーリング定数

    p(j) = floor(si(j,3) / sigma); % p = 時刻kのi
    u = si(j,3) / sigma  - p(j);

    alpha(j) = a * p(j) ^ 3 + b * p(j) ^ 2 + c * p(j) + d;

    zi(j,2) = alpha(j) - u * (alpha(j+1) - alpha(j));

    E(j+1) = E(j) + (zi(j,2) - (z3(j+1,3)- z3(j,3)) / (z1(j+1,3)- z1(j,3))) ^ 2;

    if p(j+1) > p(j)

        E(j+1) = E(j+1) + alpha()
    
    end

    si(j+1,3) = si(j,3) + u2(j) * dt;
    si(j+1,1) = si(j,1) + u1(j) * cos(si(j+1,3)) * dt;
    si(j+1,2) = si(j,2) + u1(j) * sin(si(j+1,3)) * dt;

end
