function [a]=BER(s,y)
N=length(s);
a=0;
for i=1:N
    if s(i)~=y(i)
     a=a+1;
    end
end
a=a/N;