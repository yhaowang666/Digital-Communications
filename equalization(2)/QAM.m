function [s_c,s_s] = QAM(s)
N = length(s(1,:));
s_c = zeros(1,N);s_s = s_c;
for i = 1:N
    %同向分量的映射
    if s(1,i) == 0 
        if s(2,i) == 0
            s_c(i) = 1;
        else
            s_c(i) = 3;
        end
    else
        if s(2,i) == 0
            s_c(i) = -1;
        else
            s_c(i) = -3;
        end
    end
    
    %正交分量的映射
    if s(3,i) == 0 
        if s(4,i) == 0
            s_s(i) = -3;
        else
            s_s(i) = -1;
        end
    else
        if s(4,i) == 0
            s_s(i) = 3;
        else
            s_s(i) = 1;
        end        
    end
end
    