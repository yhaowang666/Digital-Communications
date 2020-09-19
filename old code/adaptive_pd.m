function [r_c,r_s]=adaptive_pd(s_c,s_s,H,l)
N=length(l);
d=zeros(1,N);     %%d代表各个最小欧氏距离比值,也是噪声功率谱密度的比值
D=0;
for m=1:N          %%求各个星座点的最小欧式距离
    d(m)=2*abs(H(m));
    for k=1:N
        if k~=m&&d(m)>abs(H(m)-H(k))
            d(m)=abs(H(m)-H(k));
        end
    end
    D=D+d(m);
end
r=ones(1,N);
p=r./d;    %功率之比
P=0;
for m=1:N
    P=P+p(m);
end
p=p./P.*4;
n=r./p; %功率谱密度之比
%
        a1=sqrt(N0/2*n(1));
        n1=normrnd(0,a1(b),[2,N/4]);   %产生4种功率谱密度不同的双路噪声
        n_c1=n1(1,:);n_s1=n1(2,:);
        a2=sqrt(N0/2*n(2));
        n2=normrnd(0,a2(b),[2,N/4]);   
        n_c2=n(1,:);n_s2=n(2,:);
        a3=sqrt(N0/2*n(3));
        n3=normrnd(0,a3(b),[2,N/4]);   
        n_c3=n(1,:);n_s3=n(2,:);
        a4=sqrt(N0/2*n(4));
        n4=normrnd(0,a4(b),[2,N/4]);   
        n_c4=n(1,:);n_s4=n(2,:);
%
 for i=1:length(l)
         s_c(i)=s_c(i)*H(l(i));
         s_s(i)=s_s(i)*H(l(i));
         if l(i)==1
            r_c(i)=s_c(i)+n_c1(i);r_s(i)=s_s(i)+n_s1(i);
         else if l(i)==2
                 r_c(i)=s_c(i)+n_c2(i);r_s(i)=s_s(i)+n_s2(i);
             else if l(i)==3
                     r_c(i)=s_c(i)+n_c3(i);r_s(i)=s_s(i)+n_s3(i);
                 else
                     r_c(i)=s_c(i)+n_c4(i);r_s(i)=s_s(i)+n_s4(i);
                 end
             end
         end
 end
        