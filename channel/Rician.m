clc;clear all;close all;

N=1000;
s = source(N); %信源产生，序列个数为N
Eb = 1/2;
mu = 0;
SNR = 15;
N0 = Eb./(power(10,SNR/10));
sigma = sqrt(N0/2); %计算噪声的标准差

for i =1:length(sigma)
    n = normrnd(mu,sigma(i),2,N/2);   %产生服从高斯分布的双路噪声
    n_c=n(1,:);n_s=n(2,:);
    s1_c=zeros(1,N/2);s1_s=zeros(1,N/2);

    for c=1:N/2
        s1_c(c)=s(2*c-1);
        s1_s(c)=s(2*c);
    end                     %将信源分解成双路信号
    
    [s_c1,s_s1] = QPSK(s1_c,s1_s);     %进行QPSK编码
    
    B = 0.6;
    K = [5,10,15,20];
    for j = 1:length(K)
        r = normrnd(0,sqrt(1/2),2,N/2);   % 产生瑞利乘性噪声
        h = zeros(2,N/2);
        %h = ones(2,N/2).*sqrt(K(j)/(K(j)+1)) + r.*sqrt(1/(K(j)+1));
        h(1,:) = ones(1,N/2).*sqrt(K(j)/(K(j)+1)) + r(1,:).*sqrt(1/(K(j)+1));
        h(2,:) = r(2,:).*sqrt(1/(K(j)+1));
        h_i = h(1,:);h_q = h(2,:);
        s_c = s_c1.*h_i - s_s1.*h_q ;s_s = s_c1.*h_q + s_s1.*h_i;
    
        r_c = s_c + n_c;r_s = s_s + n_s;
        figure(j)
        scatter(r_c,r_s)
    	title(sprintf('Rician,SNR = %d,K = %d',SNR(i),K(j)));
    end
end

