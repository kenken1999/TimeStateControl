clear;
close all;
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
dt = 0.05; %時間刻み
Tfin = 10; %シミュレーション終了時間
xi0 = zeros(3,1); %状態ξの初期値
t1 = [0:dt:Tfin];
%u1 = ones(1,length(t1))*4; %速度
u1 = 0.1;
u2 = 0.1; %角速度
[t,xi]= ode45(@(t,xi) TwoWheel2(t,xi,t1,u1,u2),[0 Tfin],xi0,opts);
%fig = figure(1);
plot(xi(:,1),xi(:,2),"bo")
hold on;
axis equal;
grid on;