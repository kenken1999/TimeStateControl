clear;
close all;

tic



%---各格子点に関するE1の偏微分後関数 De_1 の生成----------------

s = sym('s',[2 4 2 3]); % l,m,n,iの順

l_int = sym('l_int',2); % l (k, k+1)
m_int = sym('m_int',2); % m
n_int = sym('n_int',2); % n

for j = 1 : length(k) - 1

    G = [s(2,1,1,:) - s(1,1,1,:); s(1,2,1,:) - s(1,1,1,:); s(1,1,2,:) - s(1,1,1,:)];
    H =  transpose(reshape(G,[3,3]));
    y = [si_b1(j,1) - s(1,1,1,1); si_b1(j,2) - s(1,1,1,2); si_b1(j,3) - s(1,1,1,3)];
    P = H \ y;

    for b = 1 : 3

        G2 = [s(2,b,1,:) - s(1,b,1,:); s(1,b+1,1,:) - s(1,b,1,:); s(1,b,2,:) - s(1,b,1,:)];
        H2 =  transpose(reshape(G2,[3,3]));
        y2 = [si_b1(j+1,1) - s(1,b,1,1); si_b1(j+1,2) - s(1,b,1,2); si_b1(j+1,3) - s(1,b,1,3)];
        P2 = H2 \ y2;

        e1 = ( tan(pi/12) * (m_int(1) + P(2)) - ((n_int(2) + P2(3)) - (n_int(1) + P(3))) / ((l_int(2) + P2(1)) - (l_int(1) + P(1))) ) ^ 2;

        if b == 1
            for x = 1 : 3
                De1_type1{j,1,1,1,x} = matlabFunction(diff(e1,s(1,1,1,x)), 'vars', {s(:,1:2,:,:), l_int, m_int, n_int});
                % De1_type1{j,2,1,1,x} = matlabFunction(diff(e1,s(2,1,1,x)), 'vars', {s(:,:,:,:)});
                De1_type1{j,1,2,1,x} = matlabFunction(diff(e1,s(1,2,1,x)), 'vars', {s(:,1:2,:,:), l_int, m_int, n_int});
                % De1_type1{j,1,1,2,x} = matlabFunction(diff(e1,s(1,1,2,x)), 'vars', {s(:,:,:,:)});

                %---誤差関数E計算用
                e1_type1_func{j} = matlabFunction(e1, 'vars', {s(:,1:2,:,:), l_int, m_int, n_int});
            end
        elseif b == 2
            for x = 1 : 3
                De1_type2{j,1,1,1,x} = matlabFunction(diff(e1,s(1,1,1,x)), 'vars', {s(:,1:3,:,:), l_int, m_int, n_int});
                % De1_type2{j,2,1,1,x} = matlabFunction(diff(e1,s(2,1,1,x)), 'vars', {s(:,:,:,:)});
                De1_type2{j,1,2,1,x} = matlabFunction(diff(e1,s(1,2,1,x)), 'vars', {s(:,1:3,:,:), l_int, m_int, n_int});
                % De1_type2{j,1,1,2,x} = matlabFunction(diff(e1,s(1,1,2,x)), 'vars', {s(:,:,:,:)});
                % De1_type2{j,2,2,1,x} = matlabFunction(diff(e1,s(2,2,1,x)), 'vars', {s(:,:,:,:)});
                De1_type2{j,1,3,1,x} = matlabFunction(diff(e1,s(1,3,1,x)), 'vars', {s(:,1:3,:,:), l_int, m_int, n_int});
                % De1_type2{j,1,2,2,x} = matlabFunction(diff(e1,s(1,2,2,x)), 'vars', {s(:,:,:,:)});

                %---誤差関数E計算用
                e1_type2_func{j} = matlabFunction(e1, 'vars', {s(:,1:3,:,:), l_int, m_int, n_int});          
            end     
        else
            for x = 1 : 3
                De1_type3{j,1,1,1,x} = matlabFunction(diff(e1,s(1,1,1,x)), 'vars', {s(:,1:2,:,:),s(:,3:4,:,:), l_int, m_int, n_int});
                % De1_type3{j,2,1,1,x} = matlabFunction(diff(e1,s(2,1,1,x)), 'vars', {s(:,:,:,:)});
                De1_type3{j,1,2,1,x} = matlabFunction(diff(e1,s(1,2,1,x)), 'vars', {s(:,1:2,:,:),s(:,3:4,:,:), l_int, m_int, n_int});
                % De1_type3{j,1,1,2,x} = matlabFunction(diff(e1,s(1,1,2,x)), 'vars', {s(:,:,:,:)});
                De1_type3{j,1,3,1,x} = matlabFunction(diff(e1,s(1,3,1,x)), 'vars', {s(:,1:2,:,:),s(:,3:4,:,:), l_int, m_int, n_int});
                % De1_type3{j,2,3,1,x} = matlabFunction(diff(e1,s(2,3,1,x)), 'vars', {s(:,:,:,:)});
                De1_type3{j,1,4,1,x} = matlabFunction(diff(e1,s(1,4,1,x)), 'vars', {s(:,1:2,:,:),s(:,3:4,:,:), l_int, m_int, n_int});
                % De1_type3{j,1,3,2,x} = matlabFunction(diff(e1,s(1,3,2,x)), 'vars', {s(:,:,:,:)});

                %---誤差関数E計算用
                e1_type3_func{j} = matlabFunction(e1, 'vars', {s(:,1:2,:,:),s(:,3:4,:,:), l_int, m_int, n_int});
            end     
        end

    end

end

disp("end")


%---各格子点に関するE1の偏微分後関数 De_reg の生成----------------

sr = sym('sr',[1 3 1 3]); % l,m,n,iの順

e_reg = ( ((sr(1,3,1,1) - sr(1,2,1,1)) ^ 2 + (sr(1,3,1,2) - sr(1,2,1,2)) ^ 2 + (sr(1,3,1,3) - sr(1,2,1,3)) ^ 2)...
         - ((sr(1,2,1,1) - sr(1,1,1,1)) ^ 2 + (sr(1,2,1,2) - sr(1,1,1,2)) ^ 2 + (sr(1,2,1,3) - sr(1,1,1,3)) ^ 2) ) ^ 2;

e_reg_func = matlabFunction(e_reg, 'vars', {sr(:,:,:,:)}); % 誤差関数Eの計算用

for x = 1 : 3
    De_reg{1,1,1,x} = matlabFunction(diff(e_reg,sr(1,1,1,x)), 'vars', {sr(:,:,:,:)});
    De_reg{1,2,1,x} = matlabFunction(diff(e_reg,sr(1,2,1,x)), 'vars', {sr(:,:,:,:)});
    De_reg{1,3,1,x} = matlabFunction(diff(e_reg,sr(1,3,1,x)), 'vars', {sr(:,:,:,:)});
end


% % matファイルへの保存
% save z1z2z3_estimation_func.mat


toc