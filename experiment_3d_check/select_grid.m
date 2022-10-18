function [index, index_next, index_real, rho, max_flag] = select_grid(grid_, s, t, j, index_max, index, index_next, index_real, rho, max_flag)

    break_switch = 0;
    rho_tmp = zeros(3,1);
     

    for l = 1 : index_max(1)*2-1
        for m = 1 : index_max(2)*2-1
            for n = 1 : index_max(3)*2-1

                if break_switch == 1 
                    break;
                end

                if l==index_max(1) || m==index_max(2) || n==index_max(3)
                    continue;
                end

                if l == 3
                    damy = 0;
                end
    
                if l < index_max(1)
                    l_next = l + 1;
                    index_real(j,1) = l - 1;
                    max_flag(j,1) = 1; 
                else
                    l_next = l + 1;
                    index_real(j,1) = - l + index_max(1) + 1; 
                    max_flag(j,1) = 2;
                end
    
                if m < index_max(2)
                    m_next = m + 1;
                    index_real(j,2) = m - 1;
                    max_flag(j,2) = 1;
                else
                    m_next = m + 1;
                    index_real(j,2) = - m + index_max(2) + 1;
                    max_flag(j,2) = 2;
                end
    
                if n < index_max(3)
                    n_next = n + 1;
                    index_real(j,3) = n - 1;
                    max_flag(j,3) = 1;
                else
                    n_next = n + 1;
                    index_real(j,3) = - n + index_max(3) + 1;
                    max_flag(j,3) = 2;
                end
    
                A = [grid_(l_next,m,n,:,t) - grid_(l,m,n,:,t); grid_(l,m_next,n,:,t) - grid_(l,m,n,:,t); grid_(l,m,n_next,:,t) - grid_(l,m,n,:,t)];

                B = transpose(reshape(A,[3,3]));

                x = [s(j,1) - grid_(l,m,n,1,t); s(j,2) - grid_(l,m,n,2,t); s(j,3) - grid_(l,m,n,3,t)];

                rho_tmp = B \ x;

                if (0 <= rho_tmp(1)) && (rho_tmp(1) <= 1)
                    if (0 <= rho_tmp(2)) && (rho_tmp(2) <= 1)
                        if (0 <= rho_tmp(3)) && (rho_tmp(3) <= 1)

                            rho(j,:) = rho_tmp;
    
                            index(j,1) = l;
                            index(j,2) = m;
                            index(j,3) = n;

                            index_next(j,1) = l_next;
                            index_next(j,2) = m_next;
                            index_next(j,3) = n_next;

                            break_switch = 1;
                        end
                    end
                end
            end
        end
    end

    for i = 1:3
        if max_flag(j,i) == 2
            rho(j,i) = -rho(j,i);
        end
    end

    % 欠損時の補間
    if (break_switch == 0 && j > 1)
        index(j,:) = index(j-1,:);
        index_next(j,:) = index_next(j-1,:);     
        index_real(j,:) = index_real(j-1,:);    
        rho(j,:) = rho(j-1,:);   
        disp(j)
        disp("補間しました")
    end