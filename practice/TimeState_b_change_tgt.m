clear;
close all;

dt = 0.02; %%時間刻み=離散時間Tsとして使用
Tfin = 6; %シミュレーション終了時間
t1 = [0:dt:Tfin];

u1 = ones(1,length(t1)) * (-2);
u2 = ones(1,length(t1)) * 1;

zi = zeros(length(t1),3);
zi(1,:) = [1 0 5]; %状態z = [z1, z2, z3]T の初期値を設定

v1 = u1 * cos(atan(zi(1,2)));
v2 = u2 / (cos(atan(zi(1,2))) * cos(atan(zi(1,2))));

k2 = 4;
k3 = 5;

x = zi(1,1) + 0.5 * cos(atan(zi(1,2)));
y = zi(1,3) + 0.5 * sin(atan(zi(1,2)));

hold on
axis equal;
grid on;

axis([-5 5 -3 7])

h = plot(zi(1,1),zi(1,3), 'o', 'MarkerSize' ,20, 'MarkerFaceColor', 'b');

h2 = plot(x,y,'o', 'MarkerSize' ,8, 'MarkerFaceColor', 'r');

plot(1,5,'kx','MarkerSize', 10,'LineWidth',2) % initial point
plot(2,2,'bx','MarkerSize', 10,'LineWidth',2) % changed target point
plot(0,0,'rx','MarkerSize', 10,'LineWidth',2) % origin in a state space

for i = 1:length(t1)-1

    set(h, 'XData', zi(i,1),'YData', zi(i,3));

    set(h2, 'XData', zi(i,1) + 0.5 * cos(atan(zi(i,2))),'YData', zi(i,3) + 0.5 * sin(atan(zi(i,2))));

    drawnow;

    %%%- if output gif, uncomment bellow---%%%
    % F = getframe(gcf);
    % [X,map] = rgb2ind(F.cdata,256); % RGBデータをインデックス付きデータに変更
    % if i==1
    %     % GIFファイルに書き出し
    %     imwrite(X,map,'TimeState_b_change_tgt.gif')
    %     pause(2)
    % else
    %     % 2回目以降は'append'でアニメーションを作成
    %     imwrite(X,map,'TimeState_b_change_target.gif','WriteMode','append')
    % end

    zi(i+1,1) = zi(i,1) + v1(i) * dt; %dz1 = v1(i)dt
    zi(i+1,2) = zi(i,2) + v2(i) * dt; %v2(i)/vi(1)*v1(i)*dt
    zi(i+1,3) = zi(i,3) + zi(i,2) * v1(i) * dt;

    %v1(i+1) = -0.1*zi(i+1,1); %入力v1(=u1cosθ), v1=-λz1で(λ>0の定数)z1を0に収束

    if i > 75
        % v1(i+1) = -1 * zi(i+1,1); %入力v1(=u1cosθ), v1=-λz1で(λ>0の定数)z1を0に収束
        v1(i+1) = -1 * (zi(i+1,1) - 2);
    end

    %入力v2はv1正負で場合分け
    if v1(i+1) > 0
        % v2(i+1) = -k2 * zi(i+1,2) * v1(i+1) - k3 * zi(i+1,3) * v1(i+1); 
        v2(i+1) = -k2 * zi(i+1,2) * v1(i+1) - k3 * (zi(i+1,3) - 2) * v1(i+1); 
    else
        % v2(i+1) = k2 * zi(i+1,2) * v1(i+1) - k3 * zi(i+1,3) * v1(i+1); 
        v2(i+1) = k2 * zi(i+1,2) * v1(i+1) - k3 * (zi(i+1,3) - 2) * v1(i+1); 
    end

end
