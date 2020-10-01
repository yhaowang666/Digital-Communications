clc;clear all;
close all;

N=100000;
s = source(N); %信源产生，序列个数为N
Eb = 1/2;
mu = 0;
SNR = 15;
N0 = Eb./(power(10,SNR/10));
sigma = sqrt(N0/2); %计算噪声的标准差
h1 = [0.02 0.05 0.1 -0.2 1 -0.2 0.1 0.06 0.01];
s = conv(s,h1);
for i =1:length(sigma)
    n = normrnd(mu,sigma(i),2,length(s)/2);   %产生服从高斯分布的双路噪声
    n_c=n(1,:);n_s=n(2,:);
    s1_c=zeros(1,N/2);s1_s=zeros(1,N/2);

    for c=1:length(s)/2
        s1_c(c)=s(2*c-1);
        s1_s(c)=s(2*c);
    end                     %将信源分解成双路信号
    
    [s_c1,s_s1] = QPSK(s1_c,s1_s);     %进行QPSK编码
    

    h = normrnd(0,sqrt(1/2),2,length(s)/2);              %产生瑞利乘性噪声
    h_i = h(1,:);h_q = h(2,:);
    s_c = s_c1.*h_i - s_s1.*h_q ;s_s = s_c1.*h_q + s_s1.*h_i;
    
    r_c = s_c  + n_c;r_s = s_s + n_s;
    figure(i)
    scatter(r_c,r_s)
    xlabel('In-phase');
    ylabel('Quadrature-phase');
    title(sprintf('Rayleigh,SNR = %d',SNR(i)));
end



n = 5; %抽头系数个数为2*n
c = force_zero(h1,5);
r_c1 = conv(c,r_c);r_s1 = conv(c,r_s);
figure(2)
scatter(r_c1,r_s1)

