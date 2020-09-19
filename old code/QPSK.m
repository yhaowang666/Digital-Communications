function [s_c,s_s]=QPSK(s1_c,s1_s)
N=length(s1_c);
s_c=zeros(1,N);s_s=s_c;
for i=1:N
    if s1_c(i)==0&&s1_s(i)==0
    s_c(i)=1;s_s(i)=0;
    elseif s1_c(i)==0&&s1_s(i)==1
        s_c(i)=0;s_s(i)=1;
            elseif s1_c(i)==1&&s1_s(i)==1
        s_c(i)=-1;s_s(i)=0;
            elseif s1_c(i)==1&&s1_s(i)==0
        s_c(i)=0;s_s(i)=-1;
    end
end
    