clear;
close all;

dt = 0.01;%%時間刻み=離散時間Tsとして使用
Tfin = 3;%シミュレーション終了時間
t1 = [0:dt:Tfin];

u1 = ones(1,length(t1))*4;
u2 = ones(1,length(t1))*4;

xi_z = zeros(length(t1),3);
xi_z(1,:) = [0 0 0];%状態ξの初期値を設定

h = plot(xi_z(1,1),xi_z(1,2), 'o', 'MarkerSize' ,20, 'MarkerFaceColor', 'b');

hold on;
axis equal;
grid on;
axis([-1.5 1.5 0 3])

for i = 1:length(t1)-1
    %θ(n),すなわち xi_z(i+1,3) をx(n),y(n)の計算に使う
    xi_z(i+1,3) = xi_z(i,3) + u2(i)*dt;
    xi_z(i+1,1) = xi_z(i,1) + u1(i)*cos(xi_z(i+1,3))*dt;
    xi_z(i+1,2) = xi_z(i,2) + u1(i)*sin(xi_z(i+1,3))*dt;

    set(h, 'XData', xi_z(i,1),'YData', xi_z(i,2));

    drawnow;

    F = getframe;
      % RGBデータをインデックス付きデータに変更
      [X,map] = rgb2ind(F.cdata,256);
      if i==1
          % GIFファイルに書き出し
          imwrite(X,map,'twowheels.gif')
      else
          % 2回目以降は'append'でアニメーションを作成
          imwrite(X,map,'twowheels.gif','WriteMode','append')
      end

end
