clear;
close all;

dt = 0.04;%%時間刻み=離散時間Tsとして使用
Tfin = 0.62;%シミュレーション終了時間
t1 = [0:dt:Tfin];

u1 = ones(1,length(t1)) * 5;
u2 = ones(1,length(t1)) * 5;

si = zeros(length(t1),3);
si(1,:) = [0 0 -pi/2];%状態ξの初期値を設定

x = si(1,1) + 0.15 * cos(si(1,3));
y = si(1,2) + 0.15 * sin(si(1,3));

hold on;
axis equal;
grid on;

axis([-0.5 2.5 -1.5 1.5])

h = plot(si(1,1),si(1,2), 'o', 'MarkerSize' ,20, 'MarkerFaceColor', 'b');
h2 = plot(x,y,'o', 'MarkerSize' ,8, 'MarkerFaceColor', 'r');

for i = 1:length(t1)-1

    %θ(n),すなわち si(i+1,3) をx(n),y(n)の計算に使う
    si(i+1,3) = si(i,3) + u2(i) * dt;
    si(i+1,1) = si(i,1) + u1(i) * cos(si(i+1,3)) * dt;
    si(i+1,2) = si(i,2) + u1(i) * sin(si(i+1,3)) * dt;

    set(h, 'XData', si(i,1),'YData', si(i,2));
    set(h2, 'XData', si(i,1) + 0.15 * cos(si(i,3)),'YData', si(i,2) + 0.15 * sin(si(i,3)));

    drawnow;

    F = getframe(gcf);
    % RGBデータをインデックス付きデータに変更
    [X,map] = rgb2ind(F.cdata,256);
    if i==1
        pause(2)
        % GIFファイルに書き出し
        imwrite(X,map,'twowheels.gif')
    else
        % 2回目以降は'append'でアニメーションを作成
        imwrite(X,map,'twowheels.gif','WriteMode','append')
    end

end
