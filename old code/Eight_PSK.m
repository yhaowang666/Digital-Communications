function [s_c,s_s]=Eight_PSK(s1,s2,s3)
N=length(s1);
s_c=zeros(1,N);s_s=zeros(1,N);
for i=1:N
    if s1(i)==0&&s2(i)==0&&s3(i)==0
        s_c(i)=1;s_s(i)=0;
    elseif s1(i)==0&&s2(i)==0&&s3(i)==1
        s_c(i)=sqrt(2)/2;s_s(i)=sqrt(2)/2;
    elseif s1(i)==0&&s2(i)==1&&s3(i)==1
        s_c(i)=0;s_s(i)=1;
    elseif s1(i)==0&&s2(i)==1&&s3(i)==0
        s_c(i)=-sqrt(2)/2;s_s(i)=sqrt(2)/2;
    elseif s1(i)==1&&s2(i)==1&&s3(i)==0
        s_c(i)=-1;s_s(i)=0;
    elseif s1(i)==1&&s2(i)==1&&s3(i)==1
        s_c(i)=-sqrt(2)/2;s_s(i)=-sqrt(2)/2;
    elseif s1(i)==1&&s2(i)==0&&s3(i)==1
        s_c(i)=0;s_s(i)=-1;
    else
        s_c(i)=sqrt(2)/2;s_s(i)=-sqrt(2)/2;
    end
end