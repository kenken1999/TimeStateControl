clear;
close all;

dt = 0.01; %%時間刻み=離散時間Tsとして使用
Tfin = 3; %シミュレーション終了時間
t1 = [0:dt:Tfin];

u1 = ones(1,length(t1))*4; %入力μ1, 正の定数を用いて時間z1を単調増加させる 
u2 = ones(1,length(t1))*0; %入力μ2

zi = zeros(length(t1),3);
zi(1,:) = [0 0 0]; %状態z = [z1, z2, z3]T の初期値を設定

h = plot(zi(1,1),zi(1,3), 'o', 'MarkerSize' ,20, 'MarkerFaceColor', 'b');

hold on;
axis equal;
grid on;

axis([-1.5 15 0 3])

for i = 1:length(t1)-1
    
    zi(i+1,1) = zi(i,1) + u1(i)*dt;
    zi(i+1,2) = zi(i,2) + u2(i)*u1(i)*dt;
    zi(i+1,3) = zi(i,3) + zi(i,2)*u1(i)*dt;

    set(h, 'XData', zi(i,1),'YData', zi(i,3));

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
