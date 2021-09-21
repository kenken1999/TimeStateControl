clear;
close all;

dt = 1;%%時間刻み=離散時間Tsとして使用
Tfin = 10;%シミュレーション終了時間
t1 = [0:dt:Tfin];

for i = 1 : 4

    u1_pm = (-1) ^ (i-1);

    if i < 3
        u2_pm = 1;
    else
        u2_pm = -1;
    end

    u1 = ones(length(t1),1) * u1_pm * 0.1; % 並進速度
    u2 = ones(length(t1),1) * u2_pm * (pi/32); % 回転角速度

    si = zeros(length(t1),3);
    si(1,:) = [0 0 0];%状態ξの初期値を設定

    x = si(1,1) + 0.2 * cos(si(1,3));
    y = si(1,2) + 0.2 * sin(si(1,3));

    hold on;
    axis equal;
    grid on;

    axis([-2 2 -2 2])

    h = plot(si(1,1),si(1,2), 'o', 'MarkerSize' ,20, 'MarkerFaceColor', 'b');
    h2 = plot(x,y,'o', 'MarkerSize' ,8, 'MarkerFaceColor', 'r');

    for i = 1:length(t1)

        if i ~= length(t1)
            %θ(n),すなわち si(i+1,3) をx(n),y(n)の計算に使う
            si(i+1,3) = si(i,3) + u2(i) * dt;
            si(i+1,1) = si(i,1) + u1(i) * cos(si(i+1,3)) * dt;
            si(i+1,2) = si(i,2) + u1(i) * sin(si(i+1,3)) * dt;
        end

        set(h, 'XData', si(i,1),'YData', si(i,2));
        set(h2, 'XData', si(i,1) + 0.2 * cos(si(i,3)),'YData', si(i,2) + 0.2 * sin(si(i,3)));

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

end