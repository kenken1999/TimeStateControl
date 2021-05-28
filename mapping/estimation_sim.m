clear;
close all;

dk = 0.1; %%時間刻み
Kfin = 50; %シミュレーション終了時間
k = [0:dt:Tfin];

u1 = ones(1,length(k));
u2 = ones(1,length(k));

v1 = u1 * cos(atan(zi(1,2)));
v2 = u2 / (cos(atan(zi(1,2))) * cos(atan(zi(1,2))));

si = zeros(length(t1),3); %観測するセンサ変数 , 答えは(s1, s2, s3)=(x ,y, θ)
zi = zeros(length(t1),3); %変換後の状態変数 (z1, z2, z3)=(x, tanθ, y), z1=s1, z3=s2は既知, z2=tans3は未知 

si(1,:) = [-4 5 0]; %(s1, s2, s3)=(x ,y, θ)の初期値を設定

k2 = 2;
k3 = 3;

x = zi(1,1) + 0.5 * cos(atan(zi(1,2)));
y = zi(1,3) + 0.5 * sin(atan(zi(1,2)));

hold on;
axis equal;
grid on;

axis([-5 5 -3 7])

h = plot(zi(1,1),zi(1,3), 'o', 'MarkerSize' ,15, 'MarkerFaceColor', 'b');

h2 = plot(x,y,'o', 'MarkerSize' ,5, 'MarkerFaceColor', 'r');

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

    set(h2, 'XData', zi(i,1) + 0.5 * cos(atan(zi(i,2))),'YData', zi(i,3) + 0.5 * sin(atan(zi(i,2))));

    drawnow;

    F = getframe(gcf);
      % RGBデータをインデックス付きデータに変更
      [X,map] = rgb2ind(F.cdata,256);
      if i==1
          % GIFファイルに書き出し
          imwrite(X,map,'z2_estimation.gif')
      else
          % 2回目以降は'append'でアニメーションを作成
          imwrite(X,map,'z2_estimation.gif','WriteMode','append')
      end

end
