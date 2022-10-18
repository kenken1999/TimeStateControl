function [grid_, E1_all_value, Ereg_all_value, E_all_value, e_reg_value, check_grid_m, check_DE_1, check_DE_reg] = sd_update(grid_, t, j, k, index_max, index, index_next, index_real, m_case, E1_all_value, Ereg_all_value, E_all_value, e_reg_value, De1_case1, De1_case2, De1_case3, e1_case1_func, e1_case2_func, e1_case3_func, e_reg_func, De_reg, check_grid_m, check_DE_1, check_DE_reg, s, max_flag)

    % 学習率
    eta_s1 = 0; 
    eta_s2 = 0;
    eta_s3 = 1.0 * 10 ^ (-2);

    alpha = 1;  % 正則化項E_regの割合
    
    DE_1 = zeros(index_max(2)*2, 3);  % E_1の偏微分
    DE_reg = zeros(index_max(2)*2, 3);  % E_regの偏微分

    % E_reg の値の計算
    for m = 2 : index_max(2)-1
        Ereg_all_value(t) = Ereg_all_value(t) + e_reg_func(grid_(1,m-1:m+1,1,:,t+1));

        if m == index_max(2)-1
            e_reg_value(t) = e_reg_func(grid_(1,m-1:m+1,1,:,t+1));          
        end
    end

    for m = index_max(2)+2 : index_max(2)*2-1
        Ereg_all_value(t) = Ereg_all_value(t) + e_reg_func(grid_(1,m-1:m+1,1,:,t+1));
    end

    for m = 3 : index_max(2)  % m = 1,2 は固定
        for j = 1 : length(k) - 1
            
            % 時刻kの格子点インデックスmとピッタリ一致した場合
            if m == index(j,2)
                if m_case(j) == 1
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case1{j,1,1,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case1_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));  % 誤差関数E1の作成               
                elseif m_case(j) == 2
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case2{j,1,1,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case2_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                elseif m_case(j) == 3
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case3{j,1,1,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case3_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));       
                end
            elseif m == index_next(j,2)               
                if m_case(j) == 1
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case1{j,1,2,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case1_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                elseif m_case(j) == 2
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case2{j,1,2,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case2_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                elseif m_case(j) == 3
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case3{j,1,2,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case3_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                end                        
            elseif m == index(j+1,2)         
                if m_case(j) == 3
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case3{j,1,3,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case3_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));               
                end
            elseif m == index_next(j+1,2)               
                if m_case(j) == 2
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case2{j,1,3,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case2_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                else
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case3{j,1,4,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case3_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                end
            end

            
            if j == length(k) - 1
                % DEregの生成
                if m == index_max(2)
                    for x = 1 : 3
                        DE_reg(m,x) = De_reg{1,3,1,x}(grid_(1,m-2:m,1,:,t));
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
                
                check_DE_1(t,m) = DE_1(m,3);

                % 格子点の更新
                grid_(1,m,1,1,t+1) = grid_(1,m,1,1,t) - eta_s1 * (DE_1(m,1) + alpha * DE_reg(m,1));
                grid_(1,m,1,2,t+1) = grid_(1,m,1,2,t) - eta_s2 * (DE_1(m,2) + alpha * DE_reg(m,2));

                if m == 3
                    m =  1;
                    check_grid_m(t,m) = grid_(1,m,1,3,t);
                    m =  2;
                    check_grid_m(t,m) = grid_(1,m,1,3,t);
                    m =  3;
                    check_grid_m(t,m) = grid_(1,m,1,3,t);
                else 
                    check_grid_m(t,m) = grid_(1,m,1,3,t);
                end

                grid_(1,m,1,3,t+1) = grid_(1,m,1,3,t) - eta_s3 * (DE_1(m,3) + alpha * DE_reg(m,3));
                grid_(3,m,1,3,t+1) = grid_(1,m,1,3,t+1);
                grid_(1,m,3,3,t+1) = grid_(1,m,1,3,t+1);
                grid_(3,m,3,3,t+1) = grid_(1,m,1,3,t+1);

                if grid_(1,m,1,3,t+1) <= grid_(1,m-1,1,3,t+1)
                    disp(t)
                    disp("格子の入れ替わりが発生しました");
                end

                E_all_value(t) = E1_all_value(t) + alpha * Ereg_all_value(t);
            end
        end      
    end
   

    for m = index_max(2)+3 : index_max(2)*2  % m = 1,2 は固定
        for j = 1 : length(k) - 1
            
            % 時刻kの格子点インデックスmとピッタリ一致した場合
             if m == index(j,2)
                if m_case(j) == 1
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case1{j,1,1,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case1_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));  % 誤差関数E1の作成               
                elseif m_case(j) == 2
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case2{j,1,1,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case2_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                elseif m_case(j) == 3
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case3{j,1,1,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case3_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));       
                end
            elseif m == index_next(j,2)               
                if m_case(j) == 1
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case1{j,1,2,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case1_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                elseif m_case(j) == 2
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case2{j,1,2,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case2_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                elseif m_case(j) == 3
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case3{j,1,2,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case3_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                end                        
            elseif m == index(j+1,2)         
                if m_case(j) == 3
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case3{j,1,3,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case3_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));               
                end
            elseif m == index_next(j+1,2)               
                if m_case(j) == 2
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case2{j,1,3,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case2_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                else
                    for x = 1 : 3
                        DE_1(m,x) = DE_1(m,x) + De1_case3{j,1,4,1,x,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                    end
                    E1_all_value(t) = E1_all_value(t) + e1_case3_func{j,max_flag(j,1),max_flag(j,2),max_flag(j,3)}(grid_(max_flag(j,1)*2-1:max_flag(j,1)*2,index(j,2):index_next(j,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t), grid_(1:2,index(j+1,2):index_next(j+1,2),max_flag(j,3)*2-1:max_flag(j,3)*2,:,t),index_real(j:j+1,:));
                end
            end

            if j == length(k) - 1
                % DEregの生成
                if m == index_max(2)*2
                    for x = 1 : 3
                        DE_reg(m,x) = De_reg{1,3,1,x}(grid_(1,m-2:m,1,:,t));
                    end
                elseif m == index_max(2)*2 - 1
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

                %デバック用コード
%                 check_grid_t0(t) = grid_(1,m,1,3,t);
                if m == index_max(2)+3
                    m =  index_max(2)+1;
                    check_grid_m(t,m) = grid_(1,m,1,3,t);
                    m =  index_max(2)+2;
                    check_grid_m(t,m) = grid_(1,m,1,3,t);
                    m =  index_max(2)+3;
                    check_grid_m(t,m) = grid_(1,m,1,3,t);
                else 
                    check_grid_m(t,m) = grid_(1,m,1,3,t);
                end
%                 check_DE_1(t,m) = DE_1(m,3);
%                 check_DE_reg(t,m) = DE_reg(m,3);
%                 if m < index_max(2) - 1
%                     check_DE_reg3 = [De_reg{1,1,1,x}(grid_(1,m:m+2,1,:,t)) De_reg{1,2,1,x}(grid_(1,m-1:m+1,1,:,t)) De_reg{1,3,1,x}(grid_(1,m-2:m,1,:,t))];
%                     check_grid_m = [grid_(1,m,1,1,t) grid_(1,m+1,1,1,t) grid_(1,m+2,1,1,t); grid_(1,m,1,2,t) grid_(1,m+1,1,2,t) grid_(1,m+2,1,2,t) ;grid_(1,m,1,3,t) grid_(1,m+1,1,3,t) grid_(1,m+2,1,3,t)];
%                 end

                grid_(1,m,1,3,t+1) = grid_(1,m,1,3,t) - eta_s3 * (DE_1(m,3) + alpha * DE_reg(m,3));
                grid_(3,m,1,3,t+1) = grid_(1,m,1,3,t+1);
                grid_(1,m,3,3,t+1) = grid_(1,m,1,3,t+1);
                grid_(3,m,3,3,t+1) = grid_(1,m,1,3,t+1);
                

                if grid_(1,m,1,3,t+1) >= grid_(1,m-1,1,3,t+1)
                    disp(t)
                    disp("格子の入れ替わりが発生しました");
                end

                E_all_value(t) = E1_all_value(t) + alpha * Ereg_all_value(t);
            end
        end      
    end