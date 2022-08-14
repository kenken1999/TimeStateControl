function [s,s_corr] = sampling(s, s_corr, u1, u2, k, dk)
    
    for j = 1 : length(k) - 1  
        s(j+1,3) = s(j,3) + u2(j+1) * dk;
        s(j+1,1) = s(j,1) + u1(j+1) * cos(s(j+1,3)) * dk;
        s(j+1,2) = s(j,2) + u1(j+1) * sin(s(j+1,3)) * dk;
    
        s_corr(j+1,3) = s_corr(j,3) + u2(j+1) * dk;
        s_corr(j+1,1) = s_corr(j,1) + u1(j+1) * cos(s_corr(j+1,3)) * dk;
        s_corr(j+1,2) = s_corr(j,2) + u1(j+1) * sin(s_corr(j+1,3)) * dk;    
    end