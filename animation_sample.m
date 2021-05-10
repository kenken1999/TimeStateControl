% animation_point.m

clear; close all;

% Create motion data -- (1)
t = 0:0.001:1;   % Time data
x = sin(2*pi*t); % Position data

% Draw initial figure -- (2)
figure(1);
h = plot(x(1), 0, 'o', 'MarkerSize' ,20, 'MarkerFaceColor', 'b');
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);

% Animation loop -- (3)
for i = 1:length(x)
    set(h, 'XData', x(i));
    drawnow;
    F = getframe;
      % RGBデータをインデックス付きデータに変更
      [X,map] = rgb2ind(F.cdata,256);
      if i==1
          % GIFファイルに書き出し
          imwrite(X,map,'output2.gif')
      else
          % 2回目以降は'append'でアニメーションを作成
          imwrite(X,map,'output2.gif','WriteMode','append')
      end
end
