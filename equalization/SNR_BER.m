clc;clear all;close all;
N = 1000;
s = source(N); %信源产生，序列个数为N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%16QAM

Eb = 2.5;%16QAM每个比特能量
mu = 0;
SNR = 0 : 5 : 50;
BER = zeros(1,length(SNR));
N0 = Eb./(power(10,SNR/10));
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
    
    [s_c1,s_s1] = QAM(s1);     %进行16QAM编码
    
    h1 = normrnd(0,sqrt(1/2),N/2,N/4);              %产生瑞利乘性噪声
    h_i = h1(1:N/4,:);h_q = h1(N/4+1:N/2 ,:);
    H = h_i + 1i*h_q;
    s_r = H*(s_c1 + 1i*s_s1).';

    r1 = s_r + n_c + 1i*n_s;
    
    W = inv(H'*H)*(H');
    
    W1 = W*H;
    r_ZF = W * r1;             %均衡后输出信号
    r = zeros(1,N/2);
    for j= 1:size(r_ZF,1)
        r(j) = r_ZF(j,j);
    end
    
    r_c = real(r);r_s = imag(r);   
    
    
    
    y = judgement_16QAM(r_c,r_s);     %%16QAM解码，判决输出
    BER(i) = error_rate(s,y);        %%求误比特率
end


semilogy(SNR,BER,'b*');
%axis([0 16 10^-6 0.15]);
hold on;
grid on;
xlabel('SNR/dB');ylabel('BER');
title('BER-SNR,AWGN');
% 
% BER_true = 1/4*(3*qfunc(sqrt(4/5*Eb./N0))+2*qfunc(3*sqrt(4/5*Eb./N0))+qfunc(5*sqrt(4/5*Eb./N0))); %16QAM理想误比特率
% semilogy(SNR,BER_true,'-m');
% hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QPSK
Eb = 1/2;%QPSK每个比特能量为1/2
N0 = Eb./(power(10,SNR/10));
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
    
    h1 = normrnd(0,sqrt(1/2),N,N/2);              %产生瑞利乘性噪声
    h_i = h1(1:N/2,:);h_q = h1(N/2+1:N ,:);
    H = h_i + 1i*h_q;
    s_r = H*(s_c + 1i*s_s).';

    r1 = s_r + n_c + 1i*n_s;
    
    W = inv(H'*H)*(H');
    
    W1 = W*H;
    r_ZF = W * r1;             %均衡后输出信号
    r = zeros(1,N/2);
    for j= 1:size(r_ZF,1)
        r(j) = r_ZF(j,j);
    end
    
    r_c = real(r);r_s = imag(r);
    
    y = judgement_QPSK(r_c,r_s);     %%QPSK解码，判决输出
    BER(i) = error_rate(s,y);        %%求误比特率
end


semilogy(SNR,BER,'rs');
hold on;

% BER_true = erfc(sqrt(Eb./N0))/2; %QPSK理想误比特率
% semilogy(SNR,BER_true,'-y');
% hold on

%legend('16QAM simulated','16QAM theoretical','QPSK simulated','QPSK theoretical');

