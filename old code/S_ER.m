function [a]=S_ER(s,y)
N=length(s);
a=0;
for i=1:N/3
    if s(3*i-2)~=y(3*i-2)||s(3*i-1)~=y(3*i-1)||s(3*i)~=y(3*i)
     a=a+1;
    end
end
a=a/N*3;