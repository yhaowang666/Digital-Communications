clc;clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%只对信号乘以乘性噪声的共轭，有严重的信道衰落
% N=1000;
% s = source(N); %信源产生，序列个数为N
% Eb = 1/2;
% mu = 0;
% SNR = 15;
% N0 = Eb./(power(10,SNR/10));
% sigma = sqrt(N0/2); %计算噪声的标准差
% 
% for i =1:length(sigma)
%     n = normrnd(mu,sigma(i),2,N/2);   %产生服从高斯分布的双路噪声
%     n_c=n(1,:);n_s=n(2,:);
%     s1_c=zeros(1,N/2);s1_s=zeros(1,N/2);
%     for c=1:N/2
%         s1_c(c)=s(2*c-1);
%         s1_s(c)=s(2*c);
%     end                     %将信源分解成双路信号
%     
%     [s_c1,s_s1] = QPSK(s1_c,s1_s);     %进行QPSK编码
%     
% 
%     h = normrnd(0,sqrt(1/2),2,N/2);              %产生瑞利乘性噪声
%     h_i = h(1,:);h_q = h(2,:);
%     
%     s_r = (s_c1 + 1i*s_s1).*(h_i + 1i*h_q);
%     %s_c = s_c1.*h_i - s_s1.*h_q ;s_s = s_c1.*h_q + s_s1.*h_i;
%     
%     %r_c1 = s_c  + n_c;r_s1 = s_s + n_s;
%     r_n = s_r + n_c + 1i*n_s;
%     
%     r = (h_i -1i*h_q).* r_n;
%     
%     figure(i+4)
%     %scatter(r_n)
%     scatter(real(r),imag(r))
%     xlabel('In-phase');
%     ylabel('Quadrature-phase');
%     title(sprintf('Rayleigh,SNR = %d',SNR(i)));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=1000;
s = source(N); %信源产生，序列个数为N
Eb = 1/2;
mu = 0;
SNR = 30;
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
    

    h = normrnd(0,sqrt(1/2),2,N/2);              %产生瑞利乘性噪声
    h_i = h(1,:);h_q = h(2,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     s_r = (s_c1 + 1i*s_s1).*(h_i + 1i*h_q);
% 
%     r1 = s_r + n_c + 1i*n_s;
%     
%     H = h_i + 1i*h_q;
%     W = inv(H'*H)*(H');
%     
%     W1 = W*H;
%     r_ZF = W * r1;  
%     r = zeros(1,N/2);
%     for j= 1:size(r_ZF,1)
%         r(j) = r_ZF(j,j);
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     H = h_i + 1i*h_q;     
%     W = inv(H'*H )*(H');
%     n_ZF = W * (n_c + 1i*n_s);  
%     n = zeros(1,N/2);
%     for j= 1:size(n_ZF,1)
%         n(j) = n_ZF(j,j);
%     end
%     
%     r = s_c1 + 1i*s_s1 + n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    h1 = normrnd(0,sqrt(1/2),N,N/2);              %产生瑞利乘性噪声
    h_i = h1(1:N/2,:);h_q = h1(N/2+1:N ,:);
    H = h_i + 1i*h_q;
    s_r = H*(s_c1 + 1i*s_s1).';

    r1 = s_r + n_c + 1i*n_s;
    
    W = inv(H'*H)*(H');
    
    W1 = W*H;
    r_ZF = W * r1;  
    r = zeros(1,N/2);
    for j= 1:size(r_ZF,1)
        r(j) = r_ZF(j,j);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
    scatter(real(r),imag(r))
    xlabel('In-phase');
    ylabel('Quadrature-phase');
    title(sprintf('Rayleigh,SNR = %d',SNR(i)));
end