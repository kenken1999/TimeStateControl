clear;
close all;
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

dt = 0.05; %時間刻み
Tfin = 10; %シミュレーション終了時間
xi0 = zeros(3,1); %状態ξの初期値
t1 = [0:dt:Tfin];
u1 = 0.1*t1; %速度
u2 = 0.1*t1; %角速度

for i=0:Tfin
    [t,xi]= ode45(@(t,xi) TwoWheel(t,xi,t1,u1,u2),[t1 t1+dt],xi0,opts);
    drawnow
    %fig = figure(1);
    %plot(xi(:,1),xi(:,2))
end

hold on;
axis equal;
grid on;