function [y]=judgement_QPSK(r_c,r_s)
N=length(r_c);M=2*N;
y_c=zeros(1,N);y_s=zeros(1,N);
for i=1:N
    if r_c(i)>0&&r_s(i)/r_c(i)<1&&r_s(i)/r_c(i)>-1
        y_c(i)=0;y_s(i)=0;
    else if r_s(i)>0&&(r_s(i)/r_c(i)>1||r_s(i)/r_c(i)<-1)
            y_c(i)=0;y_s(i)=1;
        else if r_c(i)<0&&r_s(i)/r_c(i)<1&&r_s(i)/r_c(i)>-1
                y_c(i)=1;y_s(i)=1;
            else
                y_c(i)=1;y_s(i)=0;
            end
        end
    end
    y(2*i-1)=y_c(i);y(2*i)=y_s(i);
end
