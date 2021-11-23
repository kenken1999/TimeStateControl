clear;
close all;

tic

%---原点付近の格子点探索・固定, およびその他初期値の決定(線形補間)------------------------------

l_max = 2;
m_max = 11;
n_max = 2;

iteration = 2000;

param_s = zeros(2*l_max-1, 2*m_max-1, 2*n_max-1, 3, iteration);

param_s(1,1,1,:,1) = [1 1 pi/4];
param_s(2,1,1,:,1) = [1+1/sqrt(2) 1+1/sqrt(2) pi/4];
param_s(1,2,1,:,1) = [1 1 pi/3];
param_s(1,1,2,:,1) = [1-1/sqrt(2) 1+1/sqrt(2) pi/4];

s_l = param_s(2,1,1,:,1) - param_s(1,1,1,:,1);
s_m = param_s(1,2,1,:,1) - param_s(1,1,1,:,1);
s_n = param_s(1,1,2,:,1) - param_s(1,1,1,:,1);


for a = 1 : 2 * l_max - 1
    for b = 1 : 2 * m_max - 1
        for c = 1 : 2 * n_max - 1

            if a <= l_max
                l_coef = a - 1;
            else
                l_coef = l_max - a;
            end

            if b <= m_max
                m_coef = b - 1;
            else
                m_coef = m_max - b;
            end

            if c <= n_max
                n_coef = c - 1;
            else
                n_coef = n_max - c;
            end

            param_s(a,b,c,:,1) = param_s(1,1,1,:,1) + l_coef * s_l + m_coef * s_m + n_coef * s_n;
    
        end
    end
end

for i = 1 : 8

    %---サンプル収集------------------------------------------------------------------------

    dk1 = 0.1;   % 時間刻み
    K1fin = 1.9;  %シミュレーション終了時間, length(k) = Kfin + 1
    k1 = [0:dk1:K1fin];

    si_b1 = zeros(length(k1),3); % 観測するセンサ変数 , s = (s1, s2, s3) = (x ,y, θ)
    si_c1 = zeros(length(k1),3); % 補正後のセンサ変数(zi,z3空間と等しい)、結果比較用

    if rem(i,8) == 1
        u1_b1 = ones(length(k1),1) * 0.5; % 並進速度
        u2_b1 = ones(length(k1),1) * 0.6; % 回転角速度
        si_b1(1,:) = [1 1 pi/4];    % (s1, s2, s3)の初期値を設定
        si_c1(1,:) = [0 0 0];
    elseif rem(i,8) == 2
        u1_b1 = ones(length(k1),1) * 0.5; % 並進速度
        u2_b1 = ones(length(k1),1) * (-0.6); % 回転角速度
        si_b1(1,:) = [1-1/sqrt(2) 1+1/sqrt(2) pi/4];    % (s1, s2, s3)の初期値を設定
        si_c1(1,:) = [0 1 0];
    elseif rem(i,8) == 3
        u1_b1 = ones(length(k1),1) * 0.5; % 並進速度
        u2_b1 = ones(length(k1),1) * 0.6; % 回転角速度
        si_b1(1,:) = [1-1/sqrt(2) 1-1/sqrt(2) pi/4];    % (s1, s2, s3)の初期値を設定
        si_c1(1,:) = [-1 0 0];
    elseif rem(i,8) == 4
        u1_b1 = ones(length(k1),1) * 0.5; % 並進速度
        u2_b1 = ones(length(k1),1) * (-0.6); % 回転角速度
        si_b1(1,:) = [1-sqrt(2) 1 pi/4];    % (s1, s2, s3)の初期値を設定
        si_c1(1,:) = [-1 1 0];
    elseif rem(i,8) == 5
        u1_b1 = ones(length(k1),1) * 0.5; % 並進速度
        u2_b1 = ones(length(k1),1) * 0.6; % 回転角速度
        si_b1(1,:) = [1 1-sqrt(2) pi/4];    % (s1, s2, s3)の初期値を設定
        si_c1(1,:) = [-1 -1 0];
    elseif rem(i,8) == 6
        u1_b1 = ones(length(k1),1) * 0.5; % 並進速度
        u2_b1 = ones(length(k1),1) * (-0.6); % 回転角速度
        si_b1(1,:) = [1-1/sqrt(2) 1-1/sqrt(2) pi/4];    % (s1, s2, s3)の初期値を設定
        si_c1(1,:) = [-1 0 0];
    elseif rem(i,8) == 7
        u1_b1 = ones(length(k1),1) * 0.5; % 並進速度
        u2_b1 = ones(length(k1),1) * 0.6; % 回転角速度
        si_b1(1,:) = [1+1/sqrt(2) 1-1/sqrt(2) pi/4];    % (s1, s2, s3)の初期値を設定
        si_c1(1,:) = [0 0 0];
    else
        u1_b1 = ones(length(k1),1) * 0.5; % 並進速度
        u2_b1 = ones(length(k1),1) * (-0.6); % 回転角速度
        si_b1(1,:) = [1 1 pi/4];    % (s1, s2, s3)の初期値を設定
        si_c1(1,:) = [0 0 0];
    end

    for j = 1 : length(k1) - 1
        
        si_b1(j+1,3) = si_b1(j,3) + u2_b1(j+1) * dk1;
        si_b1(j+1,1) = si_b1(j,1) + u1_b1(j+1) * cos(si_b1(j+1,3)) * dk1;
        si_b1(j+1,2) = si_b1(j,2) + u1_b1(j+1) * sin(si_b1(j+1,3)) * dk1;

        si_c1(j+1,3) = si_c1(j,3) + u2_b1(j+1) * dk1;
        si_c1(j+1,1) = si_c1(j,1) + u1_b1(j+1) * cos(si_c1(j+1,3)) * dk1;
        si_c1(j+1,2) = si_c1(j,2) + u1_b1(j+1) * sin(si_c1(j+1,3)) * dk1;

    end


    %---偏微分後関数の生成------------------------------------------------------------------

    %---E1の偏微分後関数の生成----------------

    s = sym('s',[2 4 2 3]); % l,m,n,iの順

    l_int = sym('l_int',2); % l (k, k+1)
    m_int = sym('m_int',2); % m
    n_int = sym('n_int',2); % n

    for j = 1 : length(k1) - 1

        G = [s(2,1,1,:) - s(1,1,1,:); s(1,2,1,:) - s(1,1,1,:); s(1,1,2,:) - s(1,1,1,:)];
        H =  transpose(reshape(G,[3,3]));
        y = [si_b1(j,1) - s(1,1,1,1); si_b1(j,2) - s(1,1,1,2); si_b1(j,3) - s(1,1,1,3)];
        P = H \ y;

        for b = 1 : 4

            if b == 4
                G2 = [s(2,b,1,:) - s(1,b,1,:); s(1,1,1,:) - s(1,b,1,:); s(1,b,2,:) - s(1,b,1,:)];
                H2 =  transpose(reshape(G2,[3,3]));
                y2 = [si_b1(j+1,1) - s(1,b,1,1); si_b1(j+1,2) - s(1,b,1,2); si_b1(j+1,3) - s(1,b,1,3)];
                P2 = H2 \ y2;
            else
                G2 = [s(2,b,1,:) - s(1,b,1,:); s(1,b+1,1,:) - s(1,b,1,:); s(1,b,2,:) - s(1,b,1,:)];
                H2 =  transpose(reshape(G2,[3,3]));
                y2 = [si_b1(j+1,1) - s(1,b,1,1); si_b1(j+1,2) - s(1,b,1,2); si_b1(j+1,3) - s(1,b,1,3)];
                P2 = H2 \ y2;
            end

            e1 = ( tan(pi/12) * (m_int(1) + P(2)) - ((n_int(2) + P2(3)) - (n_int(1) + P(3))) / ((l_int(2) + P2(1)) - (l_int(1) + P(1))) ) ^ 2;

            if b == 1
                for x = 1 : 3
                    De1_type1{i,j,1,1,1,x} = matlabFunction(diff(e1,s(1,1,1,x)), 'vars', {s(:,1:2,:,:), l_int, m_int, n_int});
                    % De1_type1{j,2,1,1,x} = matlabFunction(diff(e1,s(2,1,1,x)), 'vars', {s(:,:,:,:)});
                    De1_type1{i,j,1,2,1,x} = matlabFunction(diff(e1,s(1,2,1,x)), 'vars', {s(:,1:2,:,:), l_int, m_int, n_int});
                    % De1_type1{j,1,1,2,x} = matlabFunction(diff(e1,s(1,1,2,x)), 'vars', {s(:,:,:,:)});
                end
            elseif b == 2
                for x = 1 : 3
                    De1_type2{i,j,1,1,1,x} = matlabFunction(diff(e1,s(1,1,1,x)), 'vars', {s(:,1:3,:,:), l_int, m_int, n_int});
                    % De1_type2{j,2,1,1,x} = matlabFunction(diff(e1,s(2,1,1,x)), 'vars', {s(:,:,:,:)});
                    De1_type2{i,j,1,2,1,x} = matlabFunction(diff(e1,s(1,2,1,x)), 'vars', {s(:,1:3,:,:), l_int, m_int, n_int});
                    % De1_type2{j,1,1,2,x} = matlabFunction(diff(e1,s(1,1,2,x)), 'vars', {s(:,:,:,:)});
                    % De1_type2{j,2,2,1,x} = matlabFunction(diff(e1,s(2,2,1,x)), 'vars', {s(:,:,:,:)});
                    De1_type2{i,j,1,3,1,x} = matlabFunction(diff(e1,s(1,3,1,x)), 'vars', {s(:,1:3,:,:), l_int, m_int, n_int});
                    % De1_type2{j,1,2,2,x} = matlabFunction(diff(e1,s(1,2,2,x)), 'vars', {s(:,:,:,:)});
                end     
            elseif b == 3
                for x = 1 : 3
                    De1_type3{i,j,1,1,1,x} = matlabFunction(diff(e1,s(1,1,1,x)), 'vars', {s(:,1:2,:,:),s(:,3:4,:,:), l_int, m_int, n_int});
                    % De1_type3{j,2,1,1,x} = matlabFunction(diff(e1,s(2,1,1,x)), 'vars', {s(:,:,:,:)});
                    De1_type3{i,j,1,2,1,x} = matlabFunction(diff(e1,s(1,2,1,x)), 'vars', {s(:,1:2,:,:),s(:,3:4,:,:), l_int, m_int, n_int});
                    % De1_type3{j,1,1,2,x} = matlabFunction(diff(e1,s(1,1,2,x)), 'vars', {s(:,:,:,:)});
                    De1_type3{i,j,1,3,1,x} = matlabFunction(diff(e1,s(1,3,1,x)), 'vars', {s(:,1:2,:,:),s(:,3:4,:,:), l_int, m_int, n_int});
                    % De1_type3{j,2,3,1,x} = matlabFunction(diff(e1,s(2,3,1,x)), 'vars', {s(:,:,:,:)});
                    De1_type3{i,j,1,4,1,x} = matlabFunction(diff(e1,s(1,4,1,x)), 'vars', {s(:,1:2,:,:),s(:,3:4,:,:), l_int, m_int, n_int});
                    % De1_type3{j,1,3,2,x} = matlabFunction(diff(e1,s(1,3,2,x)), 'vars', {s(:,:,:,:)});
                end
            else
                for x = 1 : 3
                    De1_type4{i,j,1,1,1,x} = matlabFunction(diff(e1,s(1,1,1,x)), 'vars', {s(:,1:2,:,:), s(:,4,:,:), l_int, m_int, n_int});
                    % De1_type4{j,2,1,1,x} = matlabFunction(diff(e1,s(2,1,1,x)), 'vars', {s(:,:,:,:)});
                    De1_type4{i,j,1,2,1,x} = matlabFunction(diff(e1,s(1,2,1,x)), 'vars', {s(:,1:2,:,:), s(:,4,:,:), l_int, m_int, n_int});
                    % De1_type4{j,1,1,2,x} = matlabFunction(diff(e1,s(1,1,2,x)), 'vars', {s(:,:,:,:)});
                    De1_type4{i,j,1,4,1,x} = matlabFunction(diff(e1,s(1,3,1,x)), 'vars', {s(:,1:2,:,:), s(:,4,:,:), l_int, m_int, n_int});
                    % De1_type4{j,2,4,1,x} = matlabFunction(diff(e1,s(2,3,1,x)), 'vars', {s(:,:,:,:)});
                    % De1_type4{j,1,4,2,x} = matlabFunction(diff(e1,s(1,3,2,x)), 'vars', {s(:,:,:,:)});
                end     
            end

        end

    end

    disp("end")


    %---Eregの偏微分後関数の生成----------------

    sr = sym('sr',[1 3 1 3]); % l,m,n,iの順

    e_reg = ( ((sr(1,3,1,1) - sr(1,2,1,1)) ^ 2 + (sr(1,3,1,2) - sr(1,2,1,2)) ^ 2 + (sr(1,3,1,3) - sr(1,2,1,3)) ^ 2)...
            - ((sr(1,2,1,1) - sr(1,1,1,1)) ^ 2 + (sr(1,2,1,2) - sr(1,1,1,2)) ^ 2 + (sr(1,2,1,3) - sr(1,1,1,3)) ^ 2) ) ^ 2;

    for x = 1 : 3
        De_reg{1,1,1,x} = matlabFunction(diff(e_reg,sr(1,1,1,x)), 'vars', {sr(:,:,:,:)});
        De_reg{1,2,1,x} = matlabFunction(diff(e_reg,sr(1,2,1,x)), 'vars', {sr(:,:,:,:)});
        De_reg{1,3,1,x} = matlabFunction(diff(e_reg,sr(1,3,1,x)), 'vars', {sr(:,:,:,:)});
    end

end

toc