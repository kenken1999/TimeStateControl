function [grid_, E1_all_value, Ereg_all_value, E_all_value, e_reg_value] = sd_update(grid_, t, j, k, index_max, index, index_next, index_real, m_case, E1_all_value, Ereg_all_value, E_all_value, e_reg_value, De1_case1, De1_case2, De1_case3, e1_case1_func, e1_case2_func, e1_case3_func, e_reg_func, De_reg)

    % 学習率
    eta_s1 = 0; 
    eta_s2 = 0;
    eta_s3 = 1.0 * 10 ^ (-2);

    alpha = 1;  % 正則化項E_regの割合
    
    DE_1 = zeros(index_max(2), 3);  % E_1の偏微分
    DE_reg = zeros(index_max(2), 3);  % E_regの偏微分

    % E_reg の値の計算
    for m = 2 : index_max(2)-1
        Ereg_all_value(t) = Ereg_all_value(t) + e_reg_func(grid_(1,m-1:m+1,1,:,t+1));

        if m == index_max(2)-1
            e_reg_value(t) = e_reg_func(grid_(1,m-1:m+1,1,:,t+1));          
        end
    end

    for m = 3 : index_max(2)  % m = 1,2 は固定
        for j = 1 : length(k) - 1
            
            % 時刻kの格子点インデックスmとピッタリ一致した場合
            if m == index(j,2)
                if m_case(j) == 1
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case1{j,1,1,1,x}(grid_(:,index(j,2):index_next(j,2),:,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case1_func{j}(grid_(:,index(j,2):index_next(j,2),:,:,t),index_real(j:j+1,:));  % 誤差関数E1の作成               
                elseif m_case(j) == 2
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case2{j,1,1,1,x}(grid_(:,index(j,2):index_next(j+1,2),:,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case2_func{j}(grid_(:,index(j,2):index_next(j+1,2),:,:,t),index_real(j:j+1,:));
                else
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case3{j,1,1,1,x}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case3_func{j}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),index_real(j:j+1,:));       
                end
            elseif m == index_next(j,2)               
                if m_case(j) == 1
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case1{j,1,2,1,x}(grid_(:,index(j,2):index_next(j,2),:,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case1_func{j}(grid_(:,index(j,2):index_next(j,2),:,:,t),index_real(j:j+1,:));
                elseif m_case(j) == 2
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case2{j,1,2,1,x}(grid_(:,index(j,2):index_next(j+1,2),:,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case2_func{j}(grid_(:,index(j,2):index_next(j+1,2),:,:,t),index_real(j:j+1,:));
                else
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case3{j,1,2,1,x}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case3_func{j}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),index_real(j:j+1,:));
                end                         
            elseif m == index(j+1,2)               
                for x = 1 : 3
                    DE_1(m,x) = DE_1(m,x) + De1_case3{j,1,3,1,x}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),index_real(j:j+1,:));
                end
                E1_all_value(t) = E1_all_value(t) + e1_case3_func{j}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),index_real(j:j+1,:));               
            elseif m == index_next(j+1,2)               
                if m_case(j) == 2
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case2{j,1,3,1,x}(grid_(:,index(j,2):index_next(j+1,2),:,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case2_func{j}(grid_(:,index(j,2):index_next(j+1,2),:,:,t),index_real(j:j+1,:));
                else
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case3{j,1,4,1,x}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case3_func{j}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),index_real(j:j+1,:));
                end
            end
            
            if j == length(k) - 1
                % DEregの生成
                if m == 1                    
                    for x = 1 : 3
                        DE_reg(m,x) = De_reg{1,1,1,x}(grid_(1,m:m+2,1,:,t));
                    end
                elseif m == index_max(2)
                    for x = 1 : 3
                        DE_reg(m,x) = De_reg{1,3,1,x}(grid_(1,m-2:m,1,:,t));
                    end
                elseif m == 2
                    for x = 1 : 3
                        DE_reg(m,x) = De_reg{1,1,1,x}(grid_(1,m:m+2,1,:,t)) + De_reg{1,2,1,x}(grid_(1,m-1:m+1,1,:,t));
                    end
                elseif m == index_max(2) - 1
                    for x = 1 : 3
                        DE_reg(m,x) = De_reg{1,3,1,x}(grid_(1,m-2:m,1,:,t)) + De_reg{1,2,1,x}(grid_(1,m-1:m+1,1,:,t));
                    end
                else
                    for x = 1 : 3
                        DE_reg(m,x) = De_reg{1,1,1,x}(grid_(1,m:m+2,1,:,t)) + De_reg{1,2,1,x}(grid_(1,m-1:m+1,1,:,t)) + De_reg{1,3,1,x}(grid_(1,m-2:m,1,:,t));
                    end
                end

                % 格子点の更新
                grid_(1,m,1,1,t+1) = grid_(1,m,1,1,t) - eta_s1 * (DE_1(m,1) + alpha * DE_reg(m,1));
                grid_(1,m,1,2,t+1) = grid_(1,m,1,2,t) - eta_s2 * (DE_1(m,2) + alpha * DE_reg(m,2));
                check_grid_t0 = grid_(1,m,1,3,t);
                check_DE_1 = DE_1(m,3);
                check_DE_reg = DE_reg(m,3);
                grid_(1,m,1,3,t+1) = grid_(1,m,1,3,t) - eta_s3 * (DE_1(m,3) + alpha * DE_reg(m,3));
                check_grid_t1 = grid_(1,m,1,3,t+1);

                if m == 8
                    dammy = 0;
                end
                if grid_(1,m,1,3,t+1) <= grid_(1,m-1,1,3,t+1)
                    disp(t)
                    disp("格子の入れ替わりが発生しました");
                end

                E_all_value(t) = E1_all_value(t) + alpha * Ereg_all_value(t);
            end
        end      
    end