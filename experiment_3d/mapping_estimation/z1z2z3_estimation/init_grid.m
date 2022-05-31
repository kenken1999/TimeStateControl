function grid_ = init_grid(grid_, index_max)
    
    for a = 1 : index_max(1)
        for b = 1 : index_max(2)
            for c = 1 : index_max(3)
    
                if a <= index_max(1)
                    l_coef = a - 1;
                else
                    l_coef = index_max(1) - a;
                end
    
                if b <= index_max(2)
                    m_coef = b - 1;
                else
                    m_coef = index_max(2) - b;
                end
    
                if c <= index_max(3)
                    n_coef = c - 1;
                else
                    n_coef = index_max(3) - c;
                end
    
                grid_(a,b,c,:,1) = grid_(1,1,1,:,1) + l_coef*(grid_(2,1,1,:,1)-grid_(1,1,1,:,1)) + m_coef*(grid_(1,2,1,:,1) - grid_(1,1,1,:,1)) + n_coef*(grid_(1,1,2,:,1) - grid_(1,1,1,:,1));
        
            end
        end
    end