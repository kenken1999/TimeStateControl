function sd_update()

    DE1 = zeros(index_max(2),3);  % 誤差関数E_1の偏微分   
    DEreg = zeros(index_max(2),3);  % 正則化項の偏微分

    % 誤差関数Eregの作成
    for b = 2 : index_max(2)-1
        Ereg_all_value(t) = Ereg_all_value(t) + e_reg_func(grid_(1,b-1:b+1,1,:,t+1));
    end

    for b = 3 : index_max(2)  % m = 1&2 はfix

        if b < 9
            eta_s3 = 7.5 * 10 ^ (-3);
            Ereg_coef = 1.0 * 10 ^ (-1);
        else
            eta_s3 = 7.5 * 10 ^ (-3);
            Ereg_coef = 2.0 * 10 ^ (-1);
        end

        for j = 1 : length(k) - 1
            
            % 時刻kの格子点とピッタリ一致した場合
            if b == index(j,2)

                if m_case(j) == 1

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type1{j,1,1,1,x}(grid_(:,index(j,2):index_next(j,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type1_func{j}(grid_(:,index(j,2):index_next(j,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                
                elseif m_case(j) == 2

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type2{j,1,1,1,x}(grid_(:,index(j,2):index_next(j+1,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type2_func{j}(grid_(:,index(j,2):index_next(j+1,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));

                else

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type3{j,1,1,1,x}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type3_func{j}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
            
                end

            elseif b == index_next(j,2) && b ~= index(j,2)
                
                if m_case(j) == 1

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type1{j,1,2,1,x}(grid_(:,index(j,2):index_next(j,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type1_func{j}(grid_(:,index(j,2):index_next(j,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));

                elseif m_case(j) == 2

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type2{j,1,2,1,x}(grid_(:,index(j,2):index_next(j+1,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type2_func{j}(grid_(:,index(j,2):index_next(j+1,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));

                else

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type3{j,1,2,1,x}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type3_func{j}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));

                end               
            
            elseif b == index(j+1,2) && b ~= index(j,2) && b ~= index_next(j,2)
                
                for x = 1 : 3
                    DE1(b,x) = DE1(b,x) + De1_type3{j,1,3,1,x}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                end

                % 誤差関数E1の作成
                E1_all_value(t) = E1_all_value(t) + e1_type3_func{j}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                
            elseif b == index_next(j+1,2) && b ~= index(j,2) && b ~= index_next(j,2) && b ~= index(j+1,2)
                
                if m_case(j) == 2

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type2{j,1,3,1,x}(grid_(:,index(j,2):index_next(j+1,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type2_func{j}(grid_(:,index(j,2):index_next(j+1,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));

                else

                    for x = 1 : 3
                        DE1(b,x) = DE1(b,x) + De1_type3{j,1,4,1,x}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));
                    end

                    % 誤差関数E1の作成
                    E1_all_value(t) = E1_all_value(t) + e1_type3_func{j}(grid_(:,index(j,2):index_next(j,2),:,:,t), grid_(:,index(j+1,2):index_next(j+1,2),:,:,t),l_real(j:j+1),m_real(j:j+1),n_real(j:j+1));

                end

            end
            
            if j == length(k) - 1

                % DEregの生成
                if b == 1                    
                    for x = 1 : 3
                        DEreg(b,x) = De_reg{1,1,1,x}(grid_(1,b:b+2,1,:,t));
                    end
                elseif b == index_max(2)
                    for x = 1 : 3
                        DEreg(b,x) = De_reg{1,3,1,x}(grid_(1,b-2:b,1,:,t));
                    end
                elseif b == 2
                    for x = 1 : 3
                        DEreg(b,x) = De_reg{1,1,1,x}(grid_(1,b:b+2,1,:,t)) + De_reg{1,2,1,x}(grid_(1,b-1:b+1,1,:,t));
                    end
                elseif b == index_max(2) - 1
                    for x = 1 : 3
                        DEreg(b,x) = De_reg{1,3,1,x}(grid_(1,b-2:b,1,:,t)) + De_reg{1,2,1,x}(grid_(1,b-1:b+1,1,:,t));
                    end
                else
                    for x = 1 : 3
                        DEreg(b,x) = De_reg{1,1,1,x}(grid_(1,b:b+2,1,:,t)) + De_reg{1,2,1,x}(grid_(1,b-1:b+1,1,:,t)) + De_reg{1,3,1,x}(grid_(1,b-2:b,1,:,t));
                    end
                end

                % 格子点の更新
                grid_(1,b,1,1,t+1) = grid_(1,b,1,1,t) - eta_s1 * DE1(b,1) - Ereg_coef * DEreg(b,1);
                grid_(1,b,1,2,t+1) = grid_(1,b,1,2,t) - eta_s2 * DE1(b,2) - Ereg_coef * DEreg(b,2);
                grid_(1,b,1,3,t+1) = grid_(1,b,1,3,t) - eta_s3 * DE1(b,3) - Ereg_coef * DEreg(b,3);


            end

        end
        
    end

    