theta = linspace(-pi,pi);
xc = cos(theta);
yc = -sin(theta);
plot(xc,yc);
axis equal

xt = [-1 0 1 -1];
yt = [0 0 0 0];
hold on
t = area(xt,yt); % initial flat triangle
hold off
for j = 1:length(theta)-10
    xt(2) = xc(j); % determine new vertex value
    yt(2) = yc(j);
    t.XData = xt; % update data properties
    t.YData = yt;
    drawnow% display updates
    F = getframe;
      % RGBデータをインデックス付きデータに変更
      [X,map] = rgb2ind(F.cdata,256);
      if j==1
          % GIFファイルに書き出し
          imwrite(X,map,'output.gif')
      else
          % 2回目以降は'append'でアニメーションを作成
          imwrite(X,map,'output.gif','WriteMode','append')
      end
end

