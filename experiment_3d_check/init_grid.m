function grid_ = init_grid(grid_, index_max)
    
    for l = 1 : index_max(1)*2
        for m = 1 : index_max(2)*2
            for n = 1 : index_max(3)*2
    
                if l <= index_max(1)
                    l_coef = l - 1;
                else
                    l_coef = index_max(1)+1 - l;
                end
    
                if m <= index_max(2)
                    m_coef = m - 1;
                else
                    m_coef = index_max(2)+1 - m;
                end
    
                if n <= index_max(3)
                    n_coef = n - 1;
                else
                    n_coef = index_max(3)+1 - n;
                end
    
                grid_(l,m,n,:,1) = grid_(1,1,1,:,1) + l_coef*(grid_(2,1,1,:,1)-grid_(1,1,1,:,1)) + m_coef*(grid_(1,2,1,:,1) - grid_(1,1,1,:,1)) + n_coef*(grid_(1,1,2,:,1) - grid_(1,1,1,:,1));       
            end
        end
    end
