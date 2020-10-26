clc;clear all;close all;
N = 10000000;
s = source(N); %信源产生，序列个数为N
SNR = 0 : 1 : 12;

s1 = hamming_encoding(s);  %(7,4)汉明编码
N1 = length(s1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QPSK
Eb = 1/2;%QPSK每个比特能量为1/2
N0 = Eb./(power(10,SNR/10));
sigma = sqrt(N0/2); %计算噪声的标准差
mu = 0;

for i = 1:length(sigma)
    n = normrnd(mu,sigma(i),[2,N/2]);   %产生服从高斯分布的双路噪声
    n_c = n(1,:);n_s = n(2,:);
    s1_c = zeros(1,N/2);s1_s = zeros(1,N/2);

    for c = 1:N/2
        s1_c(c) = s(2*c-1);
        s1_s(c) = s(2*c);
    end                     %将信源分解成双路信号
    
    [s_c,s_s] = QPSK(s1_c,s1_s);     %进行QPSK编码
    r_c = s_c + n_c;r_s = s_s + n_s;
    
    y = judgement_QPSK(r_c,r_s);     %%QPSK解码，判决输出
    BER(i) = error_rate(s,y);        %%求误比特率
end


semilogy(SNR,BER,'-bd');
axis([0 12 10^-6 1])
hold on;
grid on;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QPSK 硬判决
Eb = 7/8;%编码后QPSK每个比特能量为1/2，则编码前为7/8
N0 = Eb./(power(10,SNR/10));
sigma = sqrt(N0/2); %计算噪声的标准差
for i = 1:length(sigma)
    n = normrnd(mu,sigma(i),[2,N1/2]);   %产生服从高斯分布的双路噪声
    n_c = n(1,:);n_s = n(2,:);
    s1_c = zeros(1,N1/2);s1_s = zeros(1,N1/2);

    for c = 1:N1/2
        s1_c(c) = s1(2*c-1);
        s1_s(c) = s1(2*c);
    end                     %将信源分解成双路信号
    
    [s_c,s_s] = QPSK(s1_c,s1_s);     %进行QPSK编码
    r_c = s_c + n_c;r_s = s_s + n_s;
    
    y1 = judgement_QPSK(r_c,r_s);     %%QPSK解码，判决输出，硬判决，直接输出0或1
    
    y = hamming_decoding(y1);         %最小汉明距离译码
    
    BER(i) = error_rate(s,y);         %求误比特率
end

semilogy(SNR,BER,'-rs');
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QPSK 软判决
% for i = 1:length(sigma)
%     n = normrnd(mu,sigma(i),[2,N1/2]);   %产生服从高斯分布的双路噪声
%     n_c = n(1,:);n_s = n(2,:);
%     s1_c = zeros(1,N1/2);s1_s = zeros(1,N1/2);
% 
%     for c = 1:N1/2
%         s1_c(c) = s1(2*c-1);
%         s1_s(c) = s1(2*c);
%     end                     %将信源分解成双路信号
%     
%     [s_c,s_s] = QPSK(s1_c,s1_s);     %进行QPSK编码
%     r_c = s_c + n_c;r_s = s_s + n_s;
%     
%     y1 = judgement_QPSK(r_c,r_s);     %%QPSK解码，判决输出，硬判决，直接输出0或1
%     
%     y = hamming_decoding2(y1);         %最小汉明距离译码
%     
%     BER(i) = error_rate(s,y);        %%求误比特率
% end
for i = 1:length(sigma)
    n = normrnd(mu,sigma(i),[2,N1/2]);   %产生服从高斯分布的双路噪声
    n_c = n(1,:);n_s = n(2,:);
    s1_c = zeros(1,N1/2);s1_s = zeros(1,N1/2);

    for c = 1:N1/2
        s1_c(c) = s1(2*c-1);
        s1_s(c) = s1(2*c);
    end                     %将信源分解成双路信号
    
    [s_c,s_s] = QPSK(s1_c,s1_s);     %进行QPSK编码
    r_c = s_c + n_c;r_s = s_s + n_s;
    
    for j = 1 : length(r_c)
        r(j*2-1) = r_c(j);
        r(j*2) = r_s(j);
    end
    y = hamming_decoding1(r);         %最小汉明距离译码
    
    BER(i) = error_rate(s,y);        %%求误比特率
end

semilogy(SNR,BER,'-kp');
hold on;
















legend('Uncoded','Hard decision','Soft decision');

