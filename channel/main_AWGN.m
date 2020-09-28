clc;clear all;close all;

N=1000;
s = source(N); %信源产生，序列个数为N
Eb = 1/2;
mu = 0;
SNR = [5,10,15,20];
N0 = Eb./(power(10,SNR/10));
sigma = sqrt(N0/2); %计算噪声的标准差

for i =1:length(sigma)
    n = normrnd(mu,sigma(i),[2,N/2]);   %产生服从高斯分布的双路噪声
    n_c=n(1,:);n_s=n(2,:);
    s1_c=zeros(1,N/2);s1_s=zeros(1,N/2);

    for c=1:N/2
        s1_c(c)=s(2*c-1);
        s1_s(c)=s(2*c);
    end                     %将信源分解成双路信号
    
    [s_c,s_s] = QPSK(s1_c,s1_s);     %进行QPSK编码
    r_c = s_c + n_c;r_s = s_s + n_s;
    figure(i)
    scatter(r_c,r_s)
    xlabel('In-phase');
    ylabel('Quadrature-phase');
    title(sprintf('AWGN,SNR = %d',SNR(i)));
end

