function [y]=judgement_EightPSK(r_c,r_s)
N=length(r_c);M=3*N;
y1=zeros(1,N);y2=zeros(1,N);y3=y1;y=zeros(1,M);
for i=1:N
    a=(r_c(i)-1)^2+(r_s(i)-0)^2;
    b=(r_c(i)-sqrt(2)/2)^2+(r_s(i)-sqrt(2)/2)^2;
    c=(r_c(i)+0)^2+(r_s(i)-1)^2;
    d=(r_c(i)+sqrt(2)/2)^2+(r_s(i)-sqrt(2)/2)^2;
    e=(r_c(i)+1)^2+(r_s(i)-0)^2;
    f=(r_c(i)+sqrt(2)/2)^2+(r_s(i)+sqrt(2)/2)^2;
    g=(r_c(i)+0)^2+(r_s(i)+1)^2;
    h=(r_c(i)-sqrt(2)/2)^2+(r_s(i)+sqrt(2)/2)^2;
    m=min([a,b,c,d,e,f,g,h]);
    if a==m
        y1(i)=0;y2(i)=0;y3(i)=0;
    elseif b==m
            y1(i)=0;y2(i)=0;y3(i)=1;
    elseif c==m
             y1(i)=0;y2(i)=1;y3(i)=1;
    elseif d==m
             y1(i)=0;y2(i)=1;y3(i)=0;
    elseif e==m
             y1(i)=1;y2(i)=1;y3(i)=0;
    elseif f==m
             y1(i)=1;y2(i)=1;y3(i)=1;
    elseif g==m
             y1(i)=1;y2(i)=0;y3(i)=1;
    else
             y1(i)=1;y2(i)=0;y3(i)=0;

    end
    y(3*i-2)=y1(i);y(3*i-1)=y2(i);y(3*i)=y3(i);
end
