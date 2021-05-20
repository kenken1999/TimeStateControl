clear;
close all;

dt = 0.03; %%時間刻み=離散時間Tsとして使用
Tfin = 15; %シミュレーション終了時間
t1 = [0:dt:Tfin];

v1 = ones(1,length(t1))*(-1);
v2 = ones(1,length(t1));

zi = zeros(length(t1),3);
zi(1,:) = [1 0 5]; %状態z = [z1, z2, z3]T の初期値を設定

k2 = 4;
k3 = 5;

x = zi(1,1)+0.5*cos(zi(1,2));
y = zi(1,3)+0.5*sin(zi(1,2));

hold on
axis equal;
grid on;

axis([-5 5 -3 7])

h = plot(zi(1,1),zi(1,3), 'o', 'MarkerSize' ,15, 'MarkerFaceColor', 'b');

h2 = plot(x,y,'o', 'MarkerSize' ,5, 'MarkerFaceColor', 'r');

plot(1,5,'kx','MarkerSize', 10,'LineWidth',2)
plot(0,0,'rx','MarkerSize', 10,'LineWidth',2)

for i = 1:length(t1)-1

    set(h, 'XData', zi(i,1),'YData', zi(i,3));

    set(h2, 'XData', zi(i,1)+0.5*cos(zi(i,2)),'YData', zi(i,3)+0.5*sin(zi(i,2)));

    drawnow;

    F = getframe(gcf);
      % RGBデータをインデックス付きデータに変更
      [X,map] = rgb2ind(F.cdata,256);
      if i==1
          % GIFファイルに書き出し
          imwrite(X,map,'TimeState_b.gif')
      else
          % 2回目以降は'append'でアニメーションを作成
          imwrite(X,map,'TimeState_b.gif','WriteMode','append')
      end

    zi(i+1,1) = zi(i,1) + v1(i)*dt; %dz1 = v1(i)dt
    zi(i+1,2) = zi(i,2) + v2(i)*dt; %v2(i)/vi(1)*v1(i)*dt
    zi(i+1,3) = zi(i,3) + zi(i,2)*v1(i)*dt;

    %v1(i+1) = -0.1*zi(i+1,1); %入力v1(=u1cosθ), v1=-λz1で(λ>0の定数)z1を0に収束

    if i > 75
        v1(i+1) = -0.5*zi(i+1,1); %入力v1(=u1cosθ), v1=-λz1で(λ>0の定数)z1を0に収束
    end

    %入力v2はv1正負で場合分け
    if v1(i+1) > 0
        v2(i+1) = -k2*zi(i+1,2)*v1(i+1)-k3*zi(i+1,3)*v1(i+1); 
    else
        v2(i+1) = k2*zi(i+1,2)*v1(i+1)-k3*zi(i+1,3)*v1(i+1); 
    end

end
