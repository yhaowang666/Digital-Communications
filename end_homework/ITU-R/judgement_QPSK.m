function [y]=judgement_QPSK(r_c,r_s)
N=length(r_c);M=2*N;
y_c=zeros(1,N);y_s=zeros(1,N);
for i=1:N
    a=(r_c(i)-sqrt(2)/2)^2+(r_s(i)-sqrt(2)/2)^2;
    b=(r_c(i)+sqrt(2)/2)^2+(r_s(i)-sqrt(2)/2)^2;
    c=(r_c(i)+sqrt(2)/2)^2+(r_s(i)+sqrt(2)/2)^2;
    d=(r_c(i)-sqrt(2)/2)^2+(r_s(i)+sqrt(2)/2)^2;
    e=min([a,b,c,d]);
    if a==e
        y_c(i)=0;y_s(i)=0;
    else if b==e
            y_c(i)=0;y_s(i)=1;
        else if c==e
                y_c(i)=1;y_s(i)=1;
            else
                y_c(i)=1;y_s(i)=0;
            end
        end
    end
    y(2*i-1)=y_c(i);y(2*i)=y_s(i);
end
