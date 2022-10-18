clear;
close all;
tic


%%%%% サンプル収集 %%%%%
dk = 0.1;   % サンプル刻み
K_fin = 1.8;  % サンプリング終了時間
k = [0:dk:K_fin];

u1 = ones(length(k),1) * 0.5;  % 並進速度
u2 = ones(length(k),1) * (0.6);  % 回転角速度

s = zeros(length(k),3);  % センサ変数 s = (s1,s2,s3) = (x,x,θ)
% s(1,:) = [1 1 pi/4];     % 初期観測
s(1,:) = [1 1+exp(-30) pi/4+exp(-30)];     %% 初期観測 %%%%%%%%%%%%%%

s_corr = zeros(length(k),3);  % 補正後のセンサ変数(z1,z3空間と等しい)、結果比較用
s_corr(1,:) = [0 0 0];

[s,s_corr] = sampling(s, s_corr, u1, u2, k, dk);  % sampling.m 呼び出し


%%%%% 各格子点に関する E_1 の偏微分 De_1 の場合分け %%%%%
sigma = [1 tan(15*pi/180) 1];  % スケーリング係数

p_sym = sym('p_sym',[2 4 2 3]);  % 格子点 p^{l,m,n}, 3次元(s空間)
index_sym = sym('index_sym', [2 3]);  % 格子点インデックス (k,k+1), {l,m,n}

for j = 1 : length(k) - 1
    A_sym = [p_sym(2,1,1,:) - p_sym(1,1,1,:); p_sym(1,2,1,:) - p_sym(1,1,1,:); p_sym(1,1,2,:) - p_sym(1,1,1,:)];
    A_sym = transpose(reshape(A_sym,[3,3]));
    x_sym = [s(j,1) - p_sym(1,1,1,1); s(j,2) - p_sym(1,1,1,2); s(j,3) - p_sym(1,1,1,3)];
    rho_sym = A_sym \ x_sym;

    for case_x = 1 : 3

        A_next = [p_sym(2,case_x,1,:) - p_sym(1,case_x,1,:); p_sym(1,case_x+1,1,:) - p_sym(1,case_x,1,:); p_sym(1,case_x,2,:) - p_sym(1,case_x,1,:)];
        A_next = transpose(reshape(A_next,[3,3]));
        x_next = [s(j+1,1) - p_sym(1,case_x,1,1); s(j+1,2) - p_sym(1,case_x,1,2); s(j+1,3) - p_sym(1,case_x,1,3)];
        rho_sym_next = A_next \ x_next;

        e1 = ( sigma(2)*(index_sym(1,2)+rho_sym(2)) - sigma(3)*((index_sym(2,3)+rho_sym_next(3))-(index_sym(1,3)+rho_sym(3))) / (sigma(1)*((index_sym(2,1)+rho_sym_next(1))-(index_sym(1,1)+rho_sym(1)))) ) ^ 2;

        if case_x == 1
            for x = 1 : 3
                De1_case1{j,1,1,1,x} = matlabFunction(diff(e1,p_sym(1,1,1,x)), 'vars', {p_sym(:,1:2,:,:), index_sym});
                De1_case1{j,1,2,1,x} = matlabFunction(diff(e1,p_sym(1,2,1,x)), 'vars', {p_sym(:,1:2,:,:), index_sym});

                e1_case1_func{j} = matlabFunction(e1, 'vars', {p_sym(:,1:2,:,:), index_sym});  % 誤差関数 E_1 計算用
            end
        elseif case_x == 2
            for x = 1 : 3
                De1_case2{j,1,1,1,x} = matlabFunction(diff(e1,p_sym(1,1,1,x)), 'vars', {p_sym(:,1:3,:,:), index_sym});
                De1_case2{j,1,2,1,x} = matlabFunction(diff(e1,p_sym(1,2,1,x)), 'vars', {p_sym(:,1:3,:,:), index_sym});
                De1_case2{j,1,3,1,x} = matlabFunction(diff(e1,p_sym(1,3,1,x)), 'vars', {p_sym(:,1:3,:,:), index_sym});

                e1_case2_func{j} = matlabFunction(e1, 'vars', {p_sym(:,1:3,:,:), index_sym});  % 誤差関数 E_1 計算用
            end     
        else
            for x = 1 : 3
                De1_case3{j,1,1,1,x} = matlabFunction(diff(e1,p_sym(1,1,1,x)), 'vars', {p_sym(:,1:2,:,:), p_sym(:,3:4,:,:), index_sym});
                De1_case3{j,1,2,1,x} = matlabFunction(diff(e1,p_sym(1,2,1,x)), 'vars', {p_sym(:,1:2,:,:), p_sym(:,3:4,:,:), index_sym});
                De1_case3{j,1,3,1,x} = matlabFunction(diff(e1,p_sym(1,3,1,x)), 'vars', {p_sym(:,1:2,:,:), p_sym(:,3:4,:,:), index_sym});
                De1_case3{j,1,4,1,x} = matlabFunction(diff(e1,p_sym(1,4,1,x)), 'vars', {p_sym(:,1:2,:,:), p_sym(:,3:4,:,:), index_sym});

                e1_case3_func{j} = matlabFunction(e1, 'vars', {p_sym(:,1:2,:,:), p_sym(:,3:4,:,:), index_sym});  % 誤差関数 E_1 計算用
            end     
        end

    end
end


%%%%% 各格子点に関する E_reg の偏微分 De_reg の場合分け %%%%%
p_sym_reg = sym('p_sym_reg',[1 3 1 3]);  % 格子点 p^{l,m,n}, 3次元

e_reg = ( sqrt((p_sym_reg(1,3,1,1) - p_sym_reg(1,2,1,1)) ^ 2 + (p_sym_reg(1,3,1,2) - p_sym_reg(1,2,1,2)) ^ 2 + (p_sym_reg(1,3,1,3) - p_sym_reg(1,2,1,3)) ^ 2) - sqrt((p_sym_reg(1,2,1,1) - p_sym_reg(1,1,1,1)) ^ 2 + (p_sym_reg(1,2,1,2) - p_sym_reg(1,1,1,2)) ^ 2 + (p_sym_reg(1,2,1,3) - p_sym_reg(1,1,1,3)) ^ 2) ) ^ 2;

e_reg_func = matlabFunction(e_reg, 'vars', {p_sym_reg(:,:,:,:)});  % 誤差関数 E_reg 計算用

for x = 1 : 3
    De_reg{1,1,1,x} = matlabFunction(diff(e_reg, p_sym_reg(1,1,1,x)), 'vars', {p_sym_reg(:,:,:,:)});
    De_reg{1,2,1,x} = matlabFunction(diff(e_reg, p_sym_reg(1,2,1,x)), 'vars', {p_sym_reg(:,:,:,:)});
    De_reg{1,3,1,x} = matlabFunction(diff(e_reg, p_sym_reg(1,3,1,x)), 'vars', {p_sym_reg(:,:,:,:)});
end


%%%%% matファイルへの保存 %%%%%
save premade_diff.mat


toc