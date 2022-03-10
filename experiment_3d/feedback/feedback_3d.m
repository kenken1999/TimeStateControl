clear;
close all;
load('../mapping_estimation/result_mat/z1z2z3_estimation.mat')
load('../mapping_estimation/result_mat/gh_estimation.mat')

%feedback_simulation------------------------------------------------------

dt = 0.01; %%時間刻み=離散時間Tsとして使用
Tfin = 3; %シミュレーション終了時間
t1 = [0:dt:Tfin];

si = zeros(length(t1),3);
si(1,:) = [1.25 2 pi/2]; %(k2, k3)=(4,5)
% si(1,:) = [1.4 1.85 5*pi/12]; %(k2, k3)=(4,5)
% si(1,:) = [1.15 2 7*pi/12]; %(k2, k3)=(5,6)

zi = zeros(length(t1),3);

u1 = ones(1,length(t1)) * (-1);
u2 = ones(1,length(t1)) * (1);

m1 = ones(1,length(t1)) * (-1); % μ1 = v1 = u1cosθ
m2 = ones(1,length(t1)); % μ2 = v2/v1 = u2/(u1*cos^3θ)

gmap = zeros(1,length(t1));
hmap = zeros(1,length(t1));

k2 = 4;
k3 = 5;

x = si(1,1) + 0.1 * cos(si(1,3));
y = si(1,2) + 0.1 * sin(si(1,3));


l_max = 2;
m_max = 11;
n_max = 2;

l_now_3 = zeros(length(t1),1);
m_now_3 = zeros(length(t1),1);
n_now_3 = zeros(length(t1),1);

l_next_3 = zeros(length(t1),1);
m_next_3 = zeros(length(t1),1);
n_next_3 = zeros(length(t1),1);

l_real_3 = zeros(length(t1),1);
m_real_3 = zeros(length(t1),1);
n_real_3 = zeros(length(t1),1);

rho_3_tmp = zeros(length(t1), l_max, m_max, n_max, 3);
rho_3 = zeros(length(t1),3);

iteration = 1224;

phi_g_now = zeros(1,imax_g);
phi_h_now = zeros(1,imax_h);



hold on;
axis equal;
grid on;

axis([0 2 0.5 2.5])

xlabel("x",'FontSize',14)
ylabel("y",'FontSize',14)

h = plot(zi(1,1),zi(1,3), 'o', 'MarkerSize' ,20, 'MarkerFaceColor', 'b');

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


for i = 1:length(t1)-1

    break_switch = 0;

    set(h, 'XData', si(i,1),'YData', si(i,2));

    set(h2, 'XData', si(i,1) + 0.1 * cos(si(i,3)),'YData', si(i,2) + 0.1 * sin(si(i,3)));

    set(gcf, 'Color', 'white'); % figureの背景を透明に設定
    set(gca, 'Color', 'white'); % axisの背景を透明に設定

    drawnow;

    xlabel("x",'FontSize',16)
    ylabel("y",'FontSize',16)

    %%%- if output gif, uncomment bellow---%%%
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

    si(i+1,3) = si(i,3) + u2(i) * dt; %観測されるs3
    si(i+1,1) = si(i,1) + u1(i) * cos(si(i+1,3)) * dt; %観測されるs1
    si(i+1,2) = si(i,2) + u1(i) * sin(si(i+1,3)) * dt; %観測されるs2


    for a = 1 : l_max

        if break_switch == 1 
            break;
        end

        for b = 1 : m_max

            if break_switch == 1 
                break;
            end

            for c = 1 : n_max

                if break_switch == 1 
                    break;
                end

                if (a == l_max) || (b == m_max) || (c == n_max)
                    continue;
                end
    
                if a < l_max
                    a2 = a + 1;
                    l_real_3(i+1) = a - 1;
                elseif a == l_max + 1
                    a2 = 1;
                    l_real_3(i+1) = -1;
                else
                    a2 = a - 1;
                    l_real_3(i+1) = - a + l_max; 
                end
    
                if b < m_max
                    b2 = b + 1;
                    m_real_3(i+1) = b - 1;
                elseif b == m_max + 1
                    b2 = 1;
                    m_real_3(i+1) = -1;
                else
                    b2 = b - 1;
                    m_real_3(i+1) = - b + m_max;
                end
    
                if c < n_max
                    c2 = c + 1;
                    n_real_3(i+1) = c - 1;
                elseif c == n_max + 1
                    c2 = 1;
                    n_real_3(i+1) = -1;
                else
                    c2 = c - 1;
                    n_real_3(i+1) = - c + n_max;
                end
    
                A = [param_s(a2,b,c,:,iteration) - param_s(a,b,c,:,iteration); param_s(a,b2,c,:,iteration) - param_s(a,b,c,:,iteration); param_s(a,b,c2,:,iteration) - param_s(a,b,c,:,iteration)];

                B = transpose(reshape(A,[3,3]));

                x = [si(i+1,1) - param_s(a,b,c,1,iteration); si(i+1,2) - param_s(a,b,c,2,iteration); si(i+1,3) - param_s(a,b,c,3,iteration)];

                rho_3_tmp(i+1,a,b,c,:) = B \ x;


                if (0 <= rho_3_tmp(i+1,a,b,c,1)) && (rho_3_tmp(i+1,a,b,c,1) <= 1)
                    if (0 <= rho_3_tmp(i+1,a,b,c,2)) && (rho_3_tmp(i+1,a,b,c,2) <= 1)
                        if (0 <= rho_3_tmp(i+1,a,b,c,3)) && (rho_3_tmp(i+1,a,b,c,3) <= 1)

                            rho_3(i+1,:) = rho_3_tmp(i+1,a,b,c,:);
    
                            l_now_3(i+1) = a;
                            m_now_3(i+1) = b;
                            n_now_3(i+1) = c;

                            l_next_3(i+1) = a2;
                            m_next_3(i+1) = b2;
                            n_next_3(i+1) = c2;

                            break_switch = 1;
    
                        end
                    end
                end

            end
        end
    end

    if break_switch == 0
        % l_now_3(i+1) = l_now_3(i);
        % m_now_3(i+1) = m_now_3(i);
        % n_now_3(i+1) = n_now_3(i);

        % l_next_3(i+1) = l_next_3(i);
        % m_next_3(i+1) = m_next_3(i);
        % n_next_3(i+1) = n_next_3(i);

        % l_real_3(i+1) = l_real_3(i);
        % m_real_3(i+1) = m_real_3(i);
        % n_real_3(i+1) = n_real_3(i);

        % a = l_now_3(i+1);
        % b = m_now_3(i+1);
        % c = n_now_3(i+1);

        % rho_3(i+1,:) = rho_3_tmp(i,a,b,c,:);
        
        disp(i)
        disp("補間しました")
    end

    zi(i+1,1) = l_real_3(i+1) + rho_3(i+1,1);
    zi(i+1,2) = tan(pi/12) * (m_real_3(i+1) + rho_3(i+1,2));
    zi(i+1,3) = n_real_3(i+1) + rho_3(i+1,3);


    for p = 1 : imax_g
        phi_g_now(p) = exp( - sqrt( (si(i+1,1) - mean_g(p,1)) ^ 2 + (si(i+1,2) - mean_g(p,2)) ^ 2 + (si(i+1,3) - mean_g(p,3)) ^ 2 ) / (2 * var_g ^ 2) );
    end
    
    for q = 1 : imax_h
        phi_h_now(q) = exp( - sqrt( (si(i+1,1) - mean_h(q,1)) ^ 2 + (si(i+1,2) - mean_h(q,2)) ^ 2 + (si(i+1,3) - mean_h(q,3)) ^ 2 ) / (2 * var_h ^ 2) );
    end

    % gmap(i+1) = ((zi(i+1,2) - zi(i,2)) * u1(i)) / ((zi(i+1,1) - zi(i,1)) * u2(i));
    % hmap(i+1) = (zi(i+1,1) - zi(i,1)) / (u1(i) * dt);

    % gmap(i+1) = 1 / (cos(si(i+1,3) - pi/4)) ^ 3;
    % hmap(i+1) = cos(si(i+1,3) - pi/4);

    gmap(i+1) = phi_g_now * w_g;
    hmap(i+1) = phi_h_now * w_h;

    if i > 1
        m1(i+1) = -2 * zi(i+1,1); %入力m1(=v1=u1cosθ), m1=-λz1で(λ>0の定数)z1を0に収束
        %m1(i+1) = -2 * (zi(i+1,1)-0.25);
    end

    %入力m2はm1正負で場合分け
    if m1(i+1) > 0
        m2(i+1) = -k2 * zi(i+1,2) - k3 * zi(i+1,3);
        % m2(i+1) = -k2 * zi(i+1,2) - k3 * (zi(i+1,3)-0.5);
    else
        m2(i+1) = k2 * zi(i+1,2) - k3 * zi(i+1,3); 
        % m2(i+1) = k2 * zi(i+1,2) - k3 * (zi(i+1,3)-0.5);
    end

    u1(i+1) = m1(i+1) / hmap(i+1);
    u2(i+1) = m2(i+1) * u1(i+1) / gmap(i+1);

end