function [y]=ML(r_c,r_s,H)
N=length(r_c);M=4*N;
y1=zeros(1,N);y2=zeros(1,N);y3=y1;y4=y1;
y=zeros(1,M);
for t=1:N
    a=(r_c(t)-H(1))^2+(r_s(t)-0)^2;
    b=(r_c(t)-H(2))^2+(r_s(t)-0)^2;
    c=(r_c(t)-H(3))^2+(r_s(t)-0)^2;
    d=(r_c(t)-H(4))^2+(r_s(t)-0)^2;
    e=(r_c(t)-0)^2+(r_s(t)-H(1))^2;
    f=(r_c(t)-0)^2+(r_s(t)-H(2))^2;
    g=(r_c(t)-0)^2+(r_s(t)-H(3))^2;
    h=(r_c(t)-0)^2+(r_s(t)-H(4))^2;
    i=(r_c(t)+H(1))^2+(r_s(t)-0)^2;
    j=(r_c(t)+H(2))^2+(r_s(t)-0)^2;
    k=(r_c(t)+H(3))^2+(r_s(t)-0)^2;
    l=(r_c(t)+H(4))^2+(r_s(t)+0)^2;
    m=(r_c(t)+0)^2+(r_s(t)+H(1))^2;
    n=(r_c(t)-0)^2+(r_s(t)+H(2))^2;
    o=(r_c(t)-0)^2+(r_s(t)+H(3))^2;
    p=(r_c(t)-0)^2+(r_s(t)+H(4))^2;
    
    w=min([a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p]);
    if a==w
        y1(t)=0;y2(t)=0;y3(t)=0;y4(t)=0;
    elseif b==w
            y1(t)=0;y2(t)=1;y3(t)=0;y4(t)=0;
    elseif c==w
             y1(t)=1;y2(t)=0;y3(t)=0;y4(t)=0;
    elseif d==w
             y1(t)=1;y2(t)=1;y3(t)=0;y4(t)=0;
    elseif e==w
             y1(t)=0;y2(t)=0;y3(t)=0;y4(t)=1;
    elseif f==w
             y1(t)=0;y2(t)=1;y3(t)=0;y4(t)=1;
    elseif g==w
             y1(t)=1;y2(t)=0;y3(t)=0;y4(t)=1;
    elseif h==w
             y1(t)=1;y2(t)=1;y3(t)=0;y4(t)=1;
    elseif i==w
             y1(t)=0;y2(t)=0;y3(t)=1;y4(t)=1;
    elseif j==w
             y1(t)=0;y2(t)=1;y3(t)=1;y4(t)=1;
    elseif k==w
             y1(t)=1;y2(t)=0;y3(t)=1;y4(t)=1;
    elseif l==w
             y1(t)=1;y2(t)=1;y3(t)=1;y4(t)=1;
    elseif m==w
             y1(t)=0;y2(t)=0;y3(t)=1;y4(t)=0;
    elseif n==w
             y1(t)=0;y2(t)=1;y3(t)=1;y4(t)=0;
    elseif o==w
             y1(t)=1;y2(t)=0;y3(t)=1;y4(t)=0;
    else
             y1(t)=1;y2(t)=1;y3(t)=1;y4(t)=0;

    end
    y(4*t-3)=y1(t);y(4*t-2)=y2(t);y(4*t-1)=y3(t);y(4*t)=y4(t);
end
