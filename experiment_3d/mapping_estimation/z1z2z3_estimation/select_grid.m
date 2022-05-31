function [index, index_next, index_real, rho, rho_tmp, break_switch] = select_grid(grid_, t, j, index_max, index, index_next, index_real, rho, rho_tmp, break_switch)

    for a = 1 : index_max(1) - 1
        for b = 1 : index_max(2) - 1
            for c = 1 : index_max(3) - 1

                if break_switch == 1 
                    break;
                end
    
                if a < index_max(1)
                    a2 = a + 1;
                    index_real(j,1) = a - 1;
                elseif a == index_max(1) + 1
                    a2 = 1;
                    index_real(j,1) = -1;
                else
                    a2 = a - 1;
                    index_real(j,1) = - a + index_max(1); 
                end
    
                if b < index_max(2)
                    b2 = b + 1;
                    index_real(j,2) = b - 1;
                elseif b == index_max(2) + 1
                    b2 = 1;
                    index_real(j,2) = -1;
                else
                    b2 = b - 1;
                    index_real(j,2) = - b + index_max(2);
                end
    
                if c < index_max(3)
                    c2 = c + 1;
                    index_real(j,3) = c - 1;
                elseif c == index_max(3) + 1
                    c2 = 1;
                    index_real(j,3) = -1;
                else
                    c2 = c - 1;
                    index_real(j,3) = - c + index_max(3);
                end
    
                A = [grid_(a2,b,c,:,t) - grid_(a,b,c,:,t); grid_(a,b2,c,:,t) - grid_(a,b,c,:,t); grid_(a,b,c2,:,t) - grid_(a,b,c,:,t)];

                B = transpose(reshape(A,[3,3]));

                x = [s(j,1) - grid_(a,b,c,1,t); s(j,2) - grid_(a,b,c,2,t); s(j,3) - grid_(a,b,c,3,t)];

                rho_tmp(j,a,b,c,:) = B \ x;

                if (0 <= rho_tmp(j,a,b,c,1)) && (rho_tmp(j,a,b,c,1) <= 1)
                    if (0 <= rho_tmp(j,a,b,c,2)) && (rho_tmp(j,a,b,c,2) <= 1)
                        if (0 <= rho_tmp(j,a,b,c,3)) && (rho_tmp(j,a,b,c,3) <= 1)

                            rho(j,:) = rho_tmp(j,a,b,c,:);
    
                            index(j,1) = a;
                            index(j,2) = b;
                            index(j,3) = c;

                            index_next(j,1) = a2;
                            index_next(j,2) = b2;
                            index_next(j,3) = c2;

                            break_switch = 1;
    
                        end
                    end
                end

            end
        end
    end