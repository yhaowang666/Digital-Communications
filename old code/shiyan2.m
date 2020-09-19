clc;clear all;close all;
a=[0,0.1,0.5,1];N=1000;e=1;M=4;
s=source(N);
for b=1:4
n=normrnd(0,a(b),[2,N/2]);
n_c=n(1,:);n_s=n(2,:);
s1_c=zeros(1,N/2);s1_s=zeros(1,N/2);
for i=1:N/2
    s1_c(i)=s(2*i-1);
    s1_s(i)=s(2*i);
end
[s_c,s_s]=QPSK(s1_c,s1_s);
r_c=s_c+n_c;r_s=s_s+n_s;
figure(b)
scatter(r_c,r_s)
grid on
end