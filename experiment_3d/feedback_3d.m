clear;
close all;
load('z_estimation.mat')
load('gh_estimation.mat')


dt = 0.01; %%時間刻み=離散時間Tsとして使用
Tfin = 3; %シミュレーション終了時間
t = [0:dt:Tfin];

s_f = zeros(length(t),3);
% s_f(1,:) = [1.25 2 pi/2]; %(k2, k3)=(3,4)
% s_f(1,:) = [1.4 1.85 5*pi/12]; %(k2, k3)=(3,4)
s_f(1,:) = [1.15 2 7*pi/12]; %(k2, k3)=(4,5)

z_f = zeros(length(t),3);

u1_f = ones(1,length(t)) * (-1);
u2_f = ones(1,length(t)) * (1);

mu1 = ones(1,length(t)) * (-1); % μ1 = v1 = u1cosθ
mu2 = ones(1,length(t)); % μ2 = v2/v1 = u2_f/(u1_f*cos^3θ)

gmap = zeros(1,length(t));
hmap = zeros(1,length(t));

k2 = 4;
k3 = 5;

x = s_f(1,1) + 0.1 * cos(s_f(1,3));
y = s_f(1,2) + 0.1 * sin(s_f(1,3));


index_f = zeros(length(t),3);
index_f_next = zeros(length(t),3);
index_f_real = zeros(length(t),3);

rho_f = zeros(length(t),3);

phi_g_now = zeros(1,imax_g);
phi_h_now = zeros(1,imax_h);

rho_f_tmp = zeros(length(t), index_max(1), index_max(2), index_max(3), 3);


hold on;
axis equal;
grid on;

axis([0 2 0.5 2.5])

xlabel("x",'FontSize',14)
ylabel("y",'FontSize',14)

%%%- if show gif, uncomment bellow---%%%
h = plot(z_f(1,1),z_f(1,3), 'o', 'MarkerSize' ,20, 'MarkerFaceColor', 'b');
h2 = plot(x,y,'o', 'MarkerSize' ,8, 'MarkerFaceColor', 'r');


plot(1+1/4,2,'kx','MarkerSize', 10,'LineWidth',2)
plot(1,1,'rx','MarkerSize', 10,'LineWidth',2)

z1_plot = 0:1:3;
z3_plot = 0:1:3;

for j = 0 : 2
    if j == 1
        plot(z1_plot, z1_plot - sqrt(2) + sqrt(2) * j, '-r')
        plot(z3_plot, - z3_plot + 2 + sqrt(2) * (j-1), '-b')
    else
        plot(z1_plot, z1_plot - sqrt(2) + sqrt(2) * j, '--k')
        plot(z3_plot, - z3_plot + 2 + sqrt(2) * (j-1), '--k')
    end
end


for i = 1:length(t)-1

    break_switch = 0;


    %%% if show the trajectory figure, uncomment bellow %%%
    % hold on;
    % plot(s_f(i,1),s_f(i,2), 'o', 'MarkerSize' ,3, 'MarkerFaceColor', 'b','MarkerEdgeColor', 'b');


    %%% if show gif, uncomment bellow %%%
    set(h, 'XData', s_f(i,1),'YData', s_f(i,2));
    set(h2, 'XData', s_f(i,1) + 0.1 * cos(s_f(i,3)),'YData', s_f(i,2) + 0.1 * sin(s_f(i,3)));
    set(gcf, 'Color', 'white'); % figureの背景を透明に設定
    set(gca, 'Color', 'white'); % axisの背景を透明に設定
    drawnow;
    xlabel("x",'FontSize',16)
    ylabel("y",'FontSize',16)


    %%% if save gif, uncomment bellow %%%
    % F = getframe(gcf);
    % % RGBデータをインデックス付きデータに変更
    % [X,map] = rgb2ind(F.cdata,256);
    % if i==1
    %     % GIFファイルに書き出し
    %     imwrite(X,map,'feedback_3dim.gif')
    % else
    %     % 2回目以降は'append'でアニメーションを作成
    %     imwrite(X,map,'feedback_3dim.gif','WriteMode','append')
    % end


    s_f(i+1,3) = s_f(i,3) + u2_f(i) * dt; %観測されるs3
    s_f(i+1,1) = s_f(i,1) + u1_f(i) * cos(s_f(i+1,3)) * dt; %観測されるs1
    s_f(i+1,2) = s_f(i,2) + u1_f(i) * sin(s_f(i+1,3)) * dt; %観測されるs2

    % 格子点の選択
    % [index_f, index_f_next, index_f_real, rho_f, break_switch] = select_grid(grid_, s_f, iteration, i, index_max, index_f, index_f_next, index_f_real, rho_f, break_switch);

    for a = 1 : index_max(1) - 1

        for b = 1 : index_max(2) - 1

            for c = 1 : index_max(3) - 1

                if break_switch == 1 
                    break;
                end

                if a < index_max(1)
                    a2 = a + 1;
                    index_f_real(i+1,1) = a - 1;
                elseif a == index_max(1) + 1
                    a2 = 1;
                    index_f_real(i+1,1) = -1;
                else
                    a2 = a - 1;
                    index_f_real(i+1,1) = - a + l_max; 
                end
    
                if b < index_max(2)
                    b2 = b + 1;
                    index_f_real(i+1,2) = b - 1;
                elseif b == index_max(2) + 1
                    b2 = 1;
                    index_f_real(i+1,2) = -1;
                else
                    b2 = b - 1;
                    index_f_real(i+1,2) = - b + m_max;
                end
    
                if c < index_max(3)
                    c2 = c + 1;
                    index_f_real(i+1,3) = c - 1;
                elseif c == index_max(3) + 1
                    c2 = 1;
                    index_f_real(i+1,3) = -1;
                else
                    c2 = c - 1;
                    index_f_real(i+1,3) = - c + n_max;
                end
    
                A = [grid_(a2,b,c,:,iteration) - grid_(a,b,c,:,iteration); grid_(a,b2,c,:,iteration) - grid_(a,b,c,:,iteration); grid_(a,b,c2,:,iteration) - grid_(a,b,c,:,iteration)];

                B = transpose(reshape(A,[3,3]));

                x = [s_f(i+1,1) - grid_(a,b,c,1,iteration); s_f(i+1,2) - grid_(a,b,c,2,iteration); s_f(i+1,3) - grid_(a,b,c,3,iteration)];

                rho_f_tmp(i+1,a,b,c,:) = B \ x;


                if (0 <= rho_f_tmp(i+1,a,b,c,1)) && (rho_f_tmp(i+1,a,b,c,1) <= 1)
                    if (0 <= rho_f_tmp(i+1,a,b,c,2)) && (rho_f_tmp(i+1,a,b,c,2) <= 1)
                        if (0 <= rho_f_tmp(i+1,a,b,c,3)) && (rho_f_tmp(i+1,a,b,c,3) <= 1)

                            rho_f(i+1,:) = rho_f_tmp(i+1,a,b,c,:);
    
                            index_f(i+1,1) = a;
                            index_f(i+1,2) = b;
                            index_f(i+1,3) = c;

                            index_f_next(i+1,1) = a2;
                            index_f_next(i+1,2) = b2;
                            index_f_next(i+1,3) = c2;

                            break_switch = 1;
    
                        end
                    end
                end

            end
        end
    end


    if break_switch == 0
        disp(i)
        disp("補間しました")
    end

    z_f(i+1,1) = sigma(1) * (index_f_real(i+1,1) + rho_f(i+1,1));
    z_f(i+1,2) = sigma(2) * (index_f_real(i+1,2) + rho_f(i+1,2));
    z_f(i+1,3) = sigma(3) * (index_f_real(i+1,3) + rho_f(i+1,3));


    for p = 1 : imax_g
        phi_g_now(p) = exp( - sqrt( (s_f(i+1,1) - mean_g(p,1)) ^ 2 + (s_f(i+1,2) - mean_g(p,2)) ^ 2 + (s_f(i+1,3) - mean_g(p,3)) ^ 2 ) / (2 * var_g ^ 2) );
    end
    
    for q = 1 : imax_h
        phi_h_now(q) = exp( - sqrt( (s_f(i+1,1) - mean_h(q,1)) ^ 2 + (s_f(i+1,2) - mean_h(q,2)) ^ 2 + (s_f(i+1,3) - mean_h(q,3)) ^ 2 ) / (2 * var_h ^ 2) );
    end

    gmap(i+1) = phi_g_now * w_g;
    hmap(i+1) = phi_h_now * w_h;

    if i > 1
        mu1(i+1) = -2 * z_f(i+1,1); %入力m1(=v1=u1cosθ), mu1=-λz1で(λ>0の定数)z1を0に収束
        %mu1(i+1) = -2 * (z_f(i+1,1)-0.25);
    end

    %入力m2はm1正負で場合分け
    if mu1(i+1) > 0
        mu2(i+1) = -k2 * z_f(i+1,2) - k3 * z_f(i+1,3);
        % mu2(i+1) = -k2 * z_f(i+1,2) - k3 * (z_f(i+1,3)-0.5);
    else
        mu2(i+1) = k2 * z_f(i+1,2) - k3 * z_f(i+1,3); 
        % mu2(i+1) = k2 * z_f(i+1,2) - k3 * (z_f(i+1,3)-0.5);
    end

    u1_f(i+1) = mu1(i+1) / hmap(i+1);
    u2_f(i+1) = mu2(i+1) * u1_f(i+1) / gmap(i+1);

end


%%% if show the feedback graph of (x,y,θ), uncomment bellow %%%

% target = ones(length(t),1);
% target2 = ones(length(t),1) * (pi/4);

% tiledlayout(3,1);

% nexttile

% plot(t, target, '--k', t, s_f(:,1), '-k','LineWidth', 2)
% ylabel("x [m]",'fontsize',16)

% nexttile

% plot(t, target, '--k', t, s_f(:,2), '-k','LineWidth', 2)
% ylabel("y [m]",'fontsize',16)

% nexttile

% plot(t, target2, '--k', t, s_f(:,3), '-k','LineWidth', 2)
% xlabel("time [s]",'fontsize',16)
% ylabel("θ [rad]",'fontsize',16)