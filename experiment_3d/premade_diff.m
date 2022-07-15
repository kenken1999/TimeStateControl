clear;
close all;

tic

%--- サンプル収集-------------------------------------------

dk = 0.1;   % サンプル刻み
K_fin = 1.8;  % サンプリング終了時間
k = [0:dk:K_fin];

u1 = ones(length(k),1) * 0.5; % 並進速度
u2 = ones(length(k),1) * (0.6); % 回転角速度

s = zeros(length(k),3); % センサ変数 s = (s1,s2,s3) = (x,y,θ)
s(1,:) = [1 1 pi/4];    % 初期観測(初期位置)

s_corr = zeros(length(k),3); % 補正後のセンサ変数(z1,z3空間と等しい)、結果比較用
s_corr(1,:) = [0 0 0];

[s,s_corr] = sampling(s, s_corr, u1, u2, k, dk); % サンプル収集


%---各格子点に関するE1の偏微分後関数 De_1 の生成----------------

grid_sym = sym('grid_sym',[2 4 2 3]); % l,m,n,iの順

l_sym = sym('l_sym',2); % l^k, l^(k+1)
m_sym = sym('m_sym',2); 
n_sym = sym('n_sym',2);

sigma = [1 tan(15*pi/180) 1];


for j = 1 : length(k) - 1

    G = [grid_sym(2,1,1,:) - grid_sym(1,1,1,:); grid_sym(1,2,1,:) - grid_sym(1,1,1,:); grid_sym(1,1,2,:) - grid_sym(1,1,1,:)];
    H = transpose(reshape(G,[3,3]));
    y = [s(j,1) - grid_sym(1,1,1,1); s(j,2) - grid_sym(1,1,1,2); s(j,3) - grid_sym(1,1,1,3)];
    rho_sym = H \ y;

    for b = 1 : 3

        G_next = [grid_sym(2,b,1,:) - grid_sym(1,b,1,:); grid_sym(1,b+1,1,:) - grid_sym(1,b,1,:); grid_sym(1,b,2,:) - grid_sym(1,b,1,:)];
        H_next = transpose(reshape(G_next,[3,3]));
        y_next = [s(j+1,1) - grid_sym(1,b,1,1); s(j+1,2) - grid_sym(1,b,1,2); s(j+1,3) - grid_sym(1,b,1,3)];
        rho_sym_next = H_next \ y_next;

        e1 = ( sigma(2)*(m_sym(1)+rho_sym(2)) - (sigma(3)*(n_sym(2)+rho_sym_next(3)) - sigma(3)*(n_sym(1)+rho_sym(3))) / (sigma(1)*(l_sym(2)+rho_sym_next(1)) - sigma(1)*(l_sym(1)+rho_sym(1))) ) ^ 2;

        if b == 1
            for x = 1 : 3
                De1_type1{j,1,1,1,x} = matlabFunction(diff(e1,grid_sym(1,1,1,x)), 'vars', {grid_sym(:,1:2,:,:), l_sym, m_sym, n_sym});
                % De1_type1{j,2,1,1,x} = matlabFunction(diff(e1,grid_sym(2,1,1,x)), 'vars', {grid_sym(:,:,:,:)});
                De1_type1{j,1,2,1,x} = matlabFunction(diff(e1,grid_sym(1,2,1,x)), 'vars', {grid_sym(:,1:2,:,:), l_sym, m_sym, n_sym});
                % De1_type1{j,1,1,2,x} = matlabFunction(diff(e1,grid_sym(1,1,2,x)), 'vars', {grid_sym(:,:,:,:)});

                % 誤差関数E計算用
                e1_type1_func{j} = matlabFunction(e1, 'vars', {grid_sym(:,1:2,:,:), l_sym, m_sym, n_sym});
            end
        elseif b == 2
            for x = 1 : 3
                De1_type2{j,1,1,1,x} = matlabFunction(diff(e1,grid_sym(1,1,1,x)), 'vars', {grid_sym(:,1:3,:,:), l_sym, m_sym, n_sym});
                % De1_type2{j,2,1,1,x} = matlabFunction(diff(e1,grid_sym(2,1,1,x)), 'vars', {grid_sym(:,:,:,:)});
                De1_type2{j,1,2,1,x} = matlabFunction(diff(e1,grid_sym(1,2,1,x)), 'vars', {grid_sym(:,1:3,:,:), l_sym, m_sym, n_sym});
                % De1_type2{j,1,1,2,x} = matlabFunction(diff(e1,grid_sym(1,1,2,x)), 'vars', {grid_sym(:,:,:,:)});
                % De1_type2{j,2,2,1,x} = matlabFunction(diff(e1,grid_sym(2,2,1,x)), 'vars', {grid_sym(:,:,:,:)});
                De1_type2{j,1,3,1,x} = matlabFunction(diff(e1,grid_sym(1,3,1,x)), 'vars', {grid_sym(:,1:3,:,:), l_sym, m_sym, n_sym});
                % De1_type2{j,1,2,2,x} = matlabFunction(diff(e1,grid_sym(1,2,2,x)), 'vars', {grid_sym(:,:,:,:)});

                % 誤差関数E計算用
                e1_type2_func{j} = matlabFunction(e1, 'vars', {grid_sym(:,1:3,:,:), l_sym, m_sym, n_sym});          
            end     
        else
            for x = 1 : 3
                De1_type3{j,1,1,1,x} = matlabFunction(diff(e1,grid_sym(1,1,1,x)), 'vars', {grid_sym(:,1:2,:,:),grid_sym(:,3:4,:,:), l_sym, m_sym, n_sym});
                % De1_type3{j,2,1,1,x} = matlabFunction(diff(e1,grid_sym(2,1,1,x)), 'vars', {grid_sym(:,:,:,:)});
                De1_type3{j,1,2,1,x} = matlabFunction(diff(e1,grid_sym(1,2,1,x)), 'vars', {grid_sym(:,1:2,:,:),grid_sym(:,3:4,:,:), l_sym, m_sym, n_sym});
                % De1_type3{j,1,1,2,x} = matlabFunction(diff(e1,grid_sym(1,1,2,x)), 'vars', {grid_sym(:,:,:,:)});
                De1_type3{j,1,3,1,x} = matlabFunction(diff(e1,grid_sym(1,3,1,x)), 'vars', {grid_sym(:,1:2,:,:),grid_sym(:,3:4,:,:), l_sym, m_sym, n_sym});
                % De1_type3{j,2,3,1,x} = matlabFunction(diff(e1,grid_sym(2,3,1,x)), 'vars', {grid_sym(:,:,:,:)});
                De1_type3{j,1,4,1,x} = matlabFunction(diff(e1,grid_sym(1,4,1,x)), 'vars', {grid_sym(:,1:2,:,:),grid_sym(:,3:4,:,:), l_sym, m_sym, n_sym});
                % De1_type3{j,1,3,2,x} = matlabFunction(diff(e1,grid_sym(1,3,2,x)), 'vars', {grid_sym(:,:,:,:)});

                % 誤差関数E計算用
                e1_type3_func{j} = matlabFunction(e1, 'vars', {grid_sym(:,1:2,:,:),grid_sym(:,3:4,:,:), l_sym, m_sym, n_sym});
            end     
        end

    end

end

disp("end")


%---各格子点に関するE1の偏微分後関数 De_reg の生成----------------

grid_reg = sym('grid_reg',[1 3 1 3]); % l,m,n,iの順

e_reg = ( sqrt((grid_reg(1,3,1,1) - grid_reg(1,2,1,1)) ^ 2 + (grid_reg(1,3,1,2) - grid_reg(1,2,1,2)) ^ 2 + (grid_reg(1,3,1,3) - grid_reg(1,2,1,3)) ^ 2) - sqrt((grid_reg(1,2,1,1) - grid_reg(1,1,1,1)) ^ 2 + (grid_reg(1,2,1,2) - grid_reg(1,1,1,2)) ^ 2 + (grid_reg(1,2,1,3) - grid_reg(1,1,1,3)) ^ 2) ) ^ 2;

e_reg_func = matlabFunction(e_reg, 'vars', {grid_reg(:,:,:,:)}); % 誤差関数Eの計算用

for x = 1 : 3
    De_reg{1,1,1,x} = matlabFunction(diff(e_reg,grid_reg(1,1,1,x)), 'vars', {grid_reg(:,:,:,:)});
    De_reg{1,2,1,x} = matlabFunction(diff(e_reg,grid_reg(1,2,1,x)), 'vars', {grid_reg(:,:,:,:)});
    De_reg{1,3,1,x} = matlabFunction(diff(e_reg,grid_reg(1,3,1,x)), 'vars', {grid_reg(:,:,:,:)});
end


% % matファイルへの保存
save premade_diff.mat


toc