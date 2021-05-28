clear;
close all;

dk = 0.1; %%時間刻み
Kfin = 50; %シミュレーション終了時間
k = [0:dt:Tfin];

u1 = ones(1,length(t1)) * 5;
u2 = ones(1,length(t1)) * 5;

si = zeros(length(t1),3); %観測するセンサ変数 , 答えは(s1, s2, s3)=(x ,y, θ)
si(1,:) = [0 0 0]; %(s1, s2, s3)=(x ,y, θ)の初期値を設定

zi = zeros(length(t1),3); %変換後の状態変数 (z1, z2, z3)=(x, tanθ, y), , z3=s2は既知, z2=tans3は未知 


for j = 1:length(k)-1

    zi(j,1) = si(j,1); %z1=s1は既知
    zi(j,3) = si(j,2); %z3=s2は既知

    si(j+1,1) = si(j,1) + u1(j) * cos(si(j+1,3)) * dt;
    si(j+1,2) = si(j,2) + u1(j) * sin(si(j+1,3)) * dt;
    si(j+1,3) = si(j,3) + u2(j) * dt;


end
