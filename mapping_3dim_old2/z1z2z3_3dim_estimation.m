clear;
close all;

tic

%---(l,m,n)=(0,0,0),(1,0,0),(0,1,0),(0,0,1)のfix, およびその他初期値の決定(線形補間)----------------------------------------------------

s = sym('s',[20 20 20 3]); % l,m,nの順

l_max = 3;
m_max = 8;
n_max = 3;

iteration = 300;

param_s = zeros(l_max, m_max, n_max, 3, iteration+1);

param_s(1,1,1,:,1) = [1 1 pi/4];
param_s(2,1,1,:,1) = [1+1/sqrt(2) 1+1/sqrt(2) pi/4];
param_s(1,2,1,:,1) = [1 1 3*pi/8];
param_s(1,1,2,:,1) = [1-1/sqrt(2) 1+1/sqrt(2) pi/4];

s_l = param_s(2,1,1,:,1) - param_s(1,1,1,:,1);
s_m = param_s(1,2,1,:,1) - param_s(1,1,1,:,1);
s_n = param_s(1,1,2,:,1) - param_s(1,1,1,:,1);


for a = 1 : l_max
    for b = 1 : m_max
        for c = 1 : n_max

            param_s(a,b,c,:,1) = param_s(1,1,1,:,1) + (a-1) * s_l + (b-1) * s_m + (c-1) * s_n;
    
        end
    end
end


%---サンプル収集と格子点の更新----------------------------------------------------

dk1 = 1;   % 時間刻み
K1fin = 15;  %シミュレーション終了時間, length(k) = Kfin + 1
k1 = [0:dk1:K1fin];

u1_b1 = ones(length(k1),1) * 0.15; % 並進速度
u2_b1 = ones(length(k1),1) * (pi/90); % 回転角速度

si_b1 = zeros(length(k1),3); % 観測するセンサ変数 , s = (s1, s2, s3) = (x ,y, θ)
si_b1(1,:) = [1 1 5*pi/12];    % (s1, s2, s3)の初期値を設定

si_c1 = zeros(length(k1),3); % 補正後のセンサ変数(zi,z3空間と等しい)、結果比較用
si_c1(1,:) = [0 0 pi/6];

l_now = zeros(length(k1),1);
m_now = zeros(length(k1),1);
n_now = zeros(length(k1),1);

rho_tmp = zeros(length(k1),3);
rho = zeros(length(k1),3);


for a = 1 : l_max - 1
    for b = 1 : m_max - 1
        for c = 1 : n_max - 1

            A = [param_s(a+1,b,c,:,1) - param_s(a,b,c,:,1); param_s(a,b+1,c,:,1) - param_s(a,b,c,:,1); param_s(a,b,c+1,:,1) - param_s(a,b,c,:,1)];

            B = transpose(reshape(A,[3,3]));

            x = [si_b1(1,1) - param_s(a,b,c,1,1); si_b1(1,2) - param_s(a,b,c,2,1); si_b1(1,3) - param_s(a,b,c,3,1)];

            rho_tmp(1,:) = B \ x;

            if (0 <= rho_tmp(1,1)) && (rho_tmp(1,1) < 1)
                if (0 <= rho_tmp(1,2)) && (rho_tmp(1,2) < 1)
                    if (0 <= rho_tmp(1,3)) && (rho_tmp(1,3) < 1)

                        rho(1,:) = rho_tmp(1,:);

                        l_now(1) = a;
                        m_now(1) = b;
                        n_now(1) = c;
                        break;

                    end
                end
            end
            
        end
    end
end

E1 = 0;


for j = 1 : length(k1) - 1
   
    si_b1(j+1,3) = si_b1(j,3) + u2_b1(j+1) * dk1;
    si_b1(j+1,1) = si_b1(j,1) + u1_b1(j+1) * cos(si_b1(j+1,3)) * dk1;
    si_b1(j+1,2) = si_b1(j,2) + u1_b1(j+1) * sin(si_b1(j+1,3)) * dk1;

    si_c1(j+1,3) = si_c1(j,3) + u2_b1(j+1) * dk1;
    si_c1(j+1,1) = si_c1(j,1) + u1_b1(j+1) * cos(si_c1(j+1,3)) * dk1;
    si_c1(j+1,2) = si_c1(j,2) + u1_b1(j+1) * sin(si_c1(j+1,3)) * dk1;


    for a = 1 : l_max - 1
        for b = 1 : m_max - 1
            for c = 1 : n_max - 1
    
                A = [param_s(a+1,b,c,:,1) - param_s(a,b,c,:,1) param_s(a,b+1,c,:,1) - param_s(a,b,c,:,1) param_s(a,b,c+1,:,1) - param_s(a,b,c,:,1)];

                B = transpose(reshape(A,[3,3]));

                x = [si_b1(j+1,1) - param_s(a,b,c,1,1); si_b1(j+1,2) - param_s(a,b,c,2,1); si_b1(j+1,3) - param_s(a,b,c,3,1)];

                rho_tmp(j+1,:) = B\x;

                if (0 <= rho_tmp(j+1,1)) && (rho_tmp(j+1,1) < 1)
                    if (0 <= rho_tmp(j+1,2)) && (rho_tmp(j+1,2) < 1)
                        if (0 <= rho_tmp(j+1,3)) && (rho_tmp(j+1,3) < 1)

                            rho(j+1,:) = rho_tmp(j+1,:);
    
                            l_now(j+1) = a;
                            m_now(j+1) = b;
                            n_now(j+1) = c;
                            break;
    
                        end
                    end
                end
                
            end
        end
    end

    l = l_now(j);
    m = m_now(j);
    n = n_now(j);

    l2 = l_now(j+1);
    m2 = l_now(j+1);
    n2 = l_now(j+1);


    G = [s(l+1,m,n,:) - s(l,m,n,:); s(l,m+1,n,:) - s(l,m,n,:); s(l,m,n+1,:) - s(l,m,n,:)];
    H =  transpose(reshape(G,[3,3]));
    y = [si_b1(j,1) - s(l,m,n,1); si_b1(j,2) - s(l,m,n,2); si_b1(j,3) - s(l,m,n,3)];
    P = H \ y;

    G2 = [s(l2+1,m2,n2,:) - s(l2,m2,n2,:); s(l2,m2+1,n2,:) - s(l2,m2,n2,:); s(l2,m2,n2+1,:) - s(l2,m2,n2,:)];
    H2 =  transpose(reshape(G2,[3,3]));
    y2 = [si_b1(j+1,1) - s(l2,m2,n2,1); si_b1(j+1,2) - s(l2,m2,n2,2); si_b1(j+1,3) - s(l2,m2,n2,3)];
    P2 = H2 \ y2;


    E1 = E1 + ( m-1 + P(2) - ((n2-1 + P2(3)) - (n-1 + P(3))) / ((l2-1 + P2(1)) - (l-1 + P(1))) ) ^ 2;

end



%---正則化項の追加-----------------------------

E4 = 0;

for a = 1 : l_max - 2
    for b = 1 : m_max - 2
        for c = 1 : n_max - 2

            E4 = E4 + (s(a+2,b,c,1) - 2 * s(a+1,b,c,1) - s(a,b,c,1)) ^ 2 + (s(a,b+2,c,1) - 2 * s(a,b+1,c,1) - s(a,b,c,1)) ^ 2 + (s(a,b,c+2,1) - 2 * s(a,b,c+1,1) - s(a,b,c,1)) ^ 2 ...
                    + (s(a+2,b,c,2) - 2 * s(a+1,b,c,2) - s(a,b,c,2)) ^ 2 + (s(a,b+2,c,2) - 2 * s(a,b+1,c,2) - s(a,b,c,2)) ^ 2 + (s(a,b,c+2,2) - 2 * s(a,b,c+1,2) - s(a,b,c,2)) ^ 2 ...
                    + (s(a+2,b,c,3) - 2 * s(a+1,b,c,3) - s(a,b,c,3)) ^ 2 + (s(a,b+2,c,3) - 2 * s(a,b+1,c,3) - s(a,b,c,3)) ^ 2 + (s(a,b,c+2,3) - 2 * s(a,b,c+1,3) - s(a,b,c,3)) ^ 2;

        end
    end
end

E1 = E1 + 0.00 * E4;



%---最急降下法によるパラメータの決定----------------------------

eta_s1 = 1.0 * 10 ^ (-8); % 学習率
eta_s2 = 1.0 * 10 ^ (-8);
eta_s3 = 1.0 * 10 ^ (-3);


%---Eの設定---------------------------

E1_initial = double(subs(E1, [s(1:l_max,1:m_max,1:n_max,:)],[param_s(1:l_max,1:m_max,1:n_max,:,1)]));
                             

disp('E1_initial = ')
disp(E1_initial)
disp('--------------------')


E1_value = zeros(1,iteration);


for a = 1:l_max
    for b = 1:m_max
        for c = 1:n_max

            DE1_s1(a,b,c) = diff(E1,s(a,b,c,1));
            DE1_s2(a,b,c) = diff(E1,s(a,b,c,2));
            DE1_s3(a,b,c) = diff(E1,s(a,b,c,3));

            DE1_s1_1{a,b,c} = matlabFunction(DE1_s1(a,b,c), 'vars', {s(1:l_max,1:m_max,1:n_max,:)});
            DE1_s2_1{a,b,c} = matlabFunction(DE1_s2(a,b,c), 'vars', {s(1:l_max,1:m_max,1:n_max,:)});
            DE1_s3_1{a,b,c} = matlabFunction(DE1_s3(a,b,c), 'vars', {s(1:l_max,1:m_max,1:n_max,:)});

            disp([a b c])
            disp('------------')

        end
    end
end



for t = 1:iteration - 1

    for a = 1:l_max
        for b = 1:m_max
            for c = 1:n_max

            DE1_s1_2 = DE1_s1_1{a,b,c}(param_s(1:l_max,1:m_max,1:n_max,:,t));
            DE1_s2_2 = DE1_s2_1{a,b,c}(param_s(1:l_max,1:m_max,1:n_max,:,t));
            DE1_s3_2 = DE1_s3_1{a,b,c}(param_s(1:l_max,1:m_max,1:n_max,:,t));

            param_s(a,b,c,1,t+1) = param_s(a,b,c,1,t) - eta_s1 * DE1_s1_2;
            param_s(a,b,c,2,t+1) = param_s(a,b,c,2,t) - eta_s2 * DE1_s2_2;
            param_s(a,b,c,3,t+1) = param_s(a,b,c,3,t) - eta_s3 * DE1_s3_2;

            end
        end
    end


    E1_value(t) = double(subs(E1, [s(1:l_max,1:m_max,1:n_max,:)],[param_s(1:l_max,1:m_max,1:n_max,:,t+1)]));


    disp('t = ')
    disp(t)
    disp('E1(t) = ')
    disp(E1_value(t))
 

    if t > 1

        if (E1_value(t) > E1_value(t-1))
            disp('Ef1_2が増加しました')
            disp('iterationを強制終了します')
            iteration = t-1;
            break;
        end
        if t == iteration - 1
            iteration = t+1;
            disp('iterationを正常に終了することができました！')
            break
        end

    end

    disp('--------------------')

end

% 推定結果の取得--------------------------------------

z1_b1 = zeros(length(k1),1);
z2_b1 = zeros(length(k1),1);
z3_b1 = zeros(length(k1),1);


for j = 1:length(k1)

    l = l_now(j);
    m = m_now(j);
    n = n_now(j);
    
    A2 = [param_s(l+1,m,n,:,iteration) - param_s(l,m,n,:,iteration) param_s(l,m+1,n,:,iteration) - param_s(l,m,n,:,iteration) param_s(l,m,n+1,:,iteration) - param_s(l,m,n,:,iteration)];
    B2 = transpose(reshape(A2,[3,3]));
    x2 = [si_b1(j,1) - param_s(l,m,n,1,iteration); si_b1(j,2) - param_s(l,m,n,2,iteration); si_b1(j,3) - param_s(l,m,n,3,iteration)];

    rho_2(j+1,:) = B2\x2;

    z1_b1(j) = l-1 + rho_2(j,1);
    z2_b1(j) = m-1 + rho_2(j,2);
    z3_b1(j) = n-1 + rho_2(j,3);

end


% 推定結果のplot--------------------------------------

figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c1(:,3), si_c1(:,1), '--m', si_c1(:,3), z1_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z1 = f1(s) = s1 の答え合わせ
xlabel("s3' = θ")
ylabel('z1 = s1')
legend('真値：s1','推定値：z1')

hold off;

figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c1(:,3), tan(si_c1(:,3)), '--m', si_c1(:,3), z2_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z1 = f1(s) = s1 の答え合わせ
xlabel("s3' = θ")
ylabel('z2 = tan(s3)')
legend('真値：tan(s3)','推定値：z2')

hold off;

figure;
hold on;
grid on;

axis([-5 5 -5 5]) % π/2 ≒ 1.57

plot(si_c1(:,3), si_c1(:,2), '--m', si_c1(:,3), z3_b1(:),'-bo','MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth', 1.5) %z3 = f3(s) = s2 の答え合わせ
xlabel("s3' = θ")
ylabel('z3 = s2')
legend('真値：s2','推定値：z3')

hold off;

toc