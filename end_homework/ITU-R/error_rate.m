function [BER]=error_rate(x,y)
N=min(length(x),length(y));
e=0;
for i=1:N
    if x(i)~=y(i)
        e=e+1;
    end
end
BER=e/N;