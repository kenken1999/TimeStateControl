clear;
close all;
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

dt = 0.05; %時間刻み
Tfin = 10; %シミュレーション終了時間
xi0 = zeros(3,1); %状態ξの初期値
t1 = [0:dt:Tfin];

for k = 0:0.05:10
    u1 = 0.1*k; %速度
    u2 = 0.1*k; %角速度
    [t,xi]= ode45(@(t,xi) TwoWheel(t,xi,t1,u1,u2),[k k+0.05],xi0,opts);
    drawnow
    %fig = figure(1);
    %plot(xi(:,1),xi(:,2))
end

hold on;
axis equal;
grid on;