clc;clear all;close all;
N = 1000000;
s = source(N); %信源产生，序列个数为N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%16QAM

Es = (2*sqrt(2)+sqrt(10));%16QAM每个符号能量
mu = 0;
SNR = -5 :1 :20;
SER = zeros(1,length(SNR));
N0 = Es./(power(10,SNR/10));
sigma = sqrt(N0/2); %计算噪声的标准差

for i = 1:length(sigma)
    n = normrnd(mu,sigma(i),[2,N/4]);   %产生服从高斯分布的双路噪声
    n_c = n(1,:);n_s = n(2,:);
    s1 = zeros(4,N/4);

    for c = 1:N/4
        s1(1,c) = s(4*c-3);
        s1(2,c) = s(4*c-2);
        s1(3,c) = s(4*c-1);
        s1(4,c) = s(4*c);
    end                     %将信源分解成四路信号
    
    [s_c,s_s] = QAM(s1);     %进行16QAM编码
    r_c = s_c + n_c;r_s = s_s + n_s;
    
    y = judgement_16QAM(r_c,r_s);     %%16QAM解码，判决输出
    SER(i) = symbol_error_16QAM(s,y);        %%求误符号率
end


semilogy(SNR,SER,'b*');
hold on;
grid on;
xlabel('SNR/dB');ylabel('BER');
title('BER-SNR,AWGN');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QPSK
Es = 1;%QPSK每个符号能量为1
N0 = Es./(power(10,SNR/10));
sigma = sqrt(N0/2); %计算噪声的标准差

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
    SER(i) = symbol_error_QPSK(s,y);        %%求误符号率
end


semilogy(SNR,SER,'rs');
hold on;
legend('16QAM','QPSK');