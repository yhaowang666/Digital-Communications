function [S]=source(L)
S=rand(1,L);
for i=1:L
    if S(i)<0.5
        S(i)=0;
    else
        S(i)=1;
    end
end
