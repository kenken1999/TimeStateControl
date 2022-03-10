clear;
close all;

dk = 0.02; %%時間刻み
Kfin = 0.62; %シミュレーション終了時間
k = [0:dk:Kfin];

u1_b = ones(1,length(k)) * 5;
u2_b = ones(1,length(k)) * 5;

si_b = zeros(length(k),3); %観測するセンサ変数 , 答えは(s1, s2, s3)=(x ,y, θ)
zi_b = zeros(length(k),3); %変換後の状態変数 (z1, z2, z3)=(x, tanθ, y), z1=s1, z3=s2は既知, z2=tans3は未知 

si_b(1,:) = [0 0 -pi/2]; %(s1, s2, s3)=(x ,y, θ)の初期値を設定

%z2_estimation----------------------------------------------------

p_now = zeros(1,length(k));

alpha = sym('alpha',[1 300]);

sigma = 0.1; %スケーリング定数
p = 1;
E = 0;

for j = 1:length(k) - 1

    si_b(j+1,3) = si_b(j,3) + u2_b(j) * dk;
    si_b(j+1,1) = si_b(j,1) + u1_b(j) * cos(si_b(j+1,3)) * dk;
    si_b(j+1,2) = si_b(j,2) + u1_b(j) * sin(si_b(j+1,3)) * dk;

    zi_b(j+1,1) = si_b(j+1,1); %z1=s1は既知
    zi_b(j+1,3) = si_b(j+1,2); %z3=s2は既知

    p_now(j+1) = floor(si_b(j+1,3) / sigma);% p = 時刻kのi, 飛び飛びor被る可能性あり
    u = si_b(j+1,3) / sigma  - p_now(j+1);

    if p_now(j+1) > p_now(j)
        p = p + 1; % 対応する格子点に番号をつけていく
    end

    %zi(j,2) = alpha(p) + u * (alpha(p+1) - alpha(p)); %z2=f(s3)

    E = E + (alpha(p) + u * (alpha(p+1) - alpha(p)) - (zi_b(j+1,3)- zi_b(j,3)) / (zi_b(j+1,1)- zi_b(j,1))) ^ 2;   

end

p_now

%正則化項の追加
% for i = 2:p

%     E = E + (alpha(i+1) - 2 * alpha(i) + alpha(i-1)) ^ 2;

% end

eta = 0.05; %学習率
iteration = 10; %パラメータ更新回数（最大）

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

    if p_now(j+1) > p_now(j)
        p = p + 1;
    end

    zi_b(j+1,2) = param(iteration,p) + u * (param(iteration,p+1) - param(iteration,p)); %z2=f(s3)

end

% hold on;
% grid on;

% axis([-1 31 -5.0 5.0])

% plot(k, tan(si(:,3)), '--', k, zi(:,2),'LineWidth', 1.5) %z2 = f(s3) = tan(s3) の答え合わせ
% xlabel('時刻 k')
% ylabel('z2 = f(s3)')
% legend('真値：tan(s3)','推定値：z2 = f(s3)')



%feedback_simulation----------------------------------------

dt = 0.01; %%時間刻み=離散時間Tsとして使用
Tfin = 6; %シミュレーション終了時間
t1 = [0:dt:Tfin];

si = zeros(length(t1),3);
si(1,:) = [1 5 0]; %(s1, s2, s3)=(x ,y, θ)の初期値を設定

zi = zeros(length(t1),3);
zi(1,:) = [1 0 5]; %(s1, s2, s3)=(x ,y, θ)の初期値を設定

u1 = ones(1,length(t1)) * (-5);
u2 = ones(1,length(t1)) * 1;

v1 = u1 * cos(atan(zi(1,2)));
v2 = u2 / (cos(atan(zi(1,2))) * cos(atan(zi(1,2))));

k2 = 4;
k3 = 5;

x = zi(1,1) + 0.5 * cos(atan(zi(1,2)));
y = zi(1,3) + 0.5 * sin(atan(zi(1,2)));

hold on;
axis equal;
grid on;

axis([-5 5 -3 7])

h = plot(zi(1,1),zi(1,3), 'o', 'MarkerSize' ,20, 'MarkerFaceColor', 'b');

h2 = plot(x,y,'o', 'MarkerSize' ,8, 'MarkerFaceColor', 'r');

plot(1,5,'kx','MarkerSize', 10,'LineWidth',2)
plot(0,0,'rx','MarkerSize', 10,'LineWidth',2)

for i = 1:length(t1)-1

    set(h, 'XData', zi(i,1),'YData', zi(i,3));

    set(h2, 'XData', zi(i,1) + 0.5 * cos(atan(zi(i,2))),'YData', zi(i,3) + 0.5 * sin(atan(zi(i,2))));

    drawnow;

    %%%- if output gif, uncomment bellow---%%%
    % F = getframe(gcf);
    %   % RGBデータをインデックス付きデータに変更
    % [X,map] = rgb2ind(F.cdata,256);
    % if i==1
    %     % GIFファイルに書き出し
    %     imwrite(X,map,'feedback_z2_estimation_b.gif')
    %     pause(2)
    % else
    %     % 2回目以降は'append'でアニメーションを作成
    %     imwrite(X,map,'feedback_z2_estimation_b.gif','WriteMode','append')
    % end

    zi(i+1,1) = zi(i,1) + v1(i) * dt; %dz1 = v1(i)dt
    %zi(i+1,2) = zi(i,2) + v2(i) * dt; %v2(i)/vi(1)*v1(i)*dt
    zi(i+1,3) = zi(i,3) + zi(i,2) * v1(i) * dt;

    si(i+1,3) = si(i,3) + v2(i) * cos(atan(zi(i,2))) * cos(atan(zi(i,2))) * dt;

    floor(si(i+1,3) / sigma)
    
    for j = 2:length(k)
        if  p_now(j) == floor(si(i+1,3) / sigma)
            zi(i+1,2) = zi_b(j,2);
        % else
        %     zi(i+1,2) = zi(i,2);
        end
    end

    %v1(i+1) = -0.1*zi(i+1,1); %入力v1(=u1cosθ), v1=-λz1で(λ>0の定数)z1を0に収束

    if i > 75
        v1(i+1) = -1 * zi(i+1,1); %入力v1(=u1cosθ), v1=-λz1で(λ>0の定数)z1を0に収束
    end

    %入力v2はv1正負で場合分け
    if v1(i+1) > 0
        v2(i+1) = -k2 * zi(i+1,2) * v1(i+1) - k3 * zi(i+1,3) * v1(i+1); 
    else
        v2(i+1) = k2 * zi(i+1,2) * v1(i+1) - k3 * zi(i+1,3) * v1(i+1); 
    end

end

zi(:,2)