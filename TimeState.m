clear;
close all;

dt = 0.1; %%時間刻み=離散時間Tsとして使用
Tfin = 50; %シミュレーション終了時間
t1 = [0:dt:Tfin];

v1 = ones(1,length(t1));
v2 = ones(1,length(t1));

zi = zeros(length(t1),3);
zi(1,:) = [-4 0 5]; %状態z = [z1, z2, z3]T の初期値を設定

k2 = 2;
k3 = 3;

h = plot(zi(1,1),zi(1,3), 'o', 'MarkerSize' ,20, 'MarkerFaceColor', 'b');

hold on;
axis equal;
grid on;

axis([-5 5 -3 7])

plot(-4,5,'kx','MarkerSize', 10,'LineWidth',2)
plot(0,0,'rx','MarkerSize', 10,'LineWidth',2)

for i = 1:length(t1)-1

    v1(i) = -0.1*zi(i,1); %入力v1(=u1cosθ), v1=-λz1で(λ>0の定数)z1を0に収束

    %入力v2はv1正負で場合分け
    if v1(i) > 0
        v2(i) = -k2*zi(i,2)*v1(i)-k3*zi(i,3)*v1(i); 
    else
        v2(i) = k2*zi(i,2)*v1(i)-k3*zi(i,3)*v1(i); 
    end
    
    zi(i+1,1) = zi(i,1) + v1(i)*dt; %dz1 = v1(i)dt
    zi(i+1,2) = zi(i,2) + v2(i)*dt; %v2(i)/vi(1)*v1(i)*dt
    zi(i+1,3) = zi(i,3) + zi(i,2)*v1(i)*dt;

    set(h, 'XData', zi(i,1),'YData', zi(i,3));

    drawnow;

    F = getframe(gcf);
      % RGBデータをインデックス付きデータに変更
      [X,map] = rgb2ind(F.cdata,256);
      if i==1
          % GIFファイルに書き出し
          imwrite(X,map,'TimeState.gif')
      else
          % 2回目以降は'append'でアニメーションを作成
          imwrite(X,map,'TimeState.gif','WriteMode','append')
      end

end
