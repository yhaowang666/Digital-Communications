clc;clear all;close all;
N = 1000000;
s = source(N); %��Դ���������и���ΪN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%16QAM

Eb = (2*sqrt(2)+sqrt(10))/4;%16QAMÿ����������
mu = 0;
SNR = 0 :1 : 16;
BER = zeros(1,length(SNR));
N0 = Eb./(power(10,SNR/10));
sigma = sqrt(N0/2); %���������ı�׼��


for i = 1:length(sigma)
    n = normrnd(mu,sigma(i),[2,N/4]);   %�������Ӹ�˹�ֲ���˫·����
    n_c = n(1,:);n_s = n(2,:);
    s1 = zeros(4,N/4);

    for c = 1:N/4
        s1(1,c) = s(4*c-3);
        s1(2,c) = s(4*c-2);
        s1(3,c) = s(4*c-1);
        s1(4,c) = s(4*c);
    end                     %����Դ�ֽ����·�ź�
    
    [s_c,s_s] = QAM(s1);     %����16QAM����
    r_c = s_c + n_c;r_s = s_s + n_s;
    
    y = judgement_16QAM(r_c,r_s);     %%16QAM���룬�о����
    BER(i) = error_rate(s,y);        %%���������
end


semilogy(SNR,BER,'b*');
axis([0 16 10^-6 0.15]);
hold on;
grid on;
xlabel('SNR/dB');ylabel('BER');
title('BER-SNR,AWGN');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QPSK
Eb = 1/2;%QPSKÿ����������Ϊ1/2
N0 = Eb./(power(10,SNR/10));
sigma = sqrt(N0/2); %���������ı�׼��

for i = 1:length(sigma)
    n = normrnd(mu,sigma(i),[2,N/2]);   %�������Ӹ�˹�ֲ���˫·����
    n_c = n(1,:);n_s = n(2,:);
    s1_c = zeros(1,N/2);s1_s = zeros(1,N/2);

    for c = 1:N/2
        s1_c(c) = s(2*c-1);
        s1_s(c) = s(2*c);
    end                     %����Դ�ֽ��˫·�ź�
    
    [s_c,s_s] = QPSK(s1_c,s1_s);     %����QPSK����
    r_c = s_c + n_c;r_s = s_s + n_s;
    
    y = judgement_QPSK(r_c,r_s);     %%QPSK���룬�о����
    BER(i) = error_rate(s,y);        %%���������
end


semilogy(SNR,BER,'rs');
hold on;

BER_true = erfc(sqrt(Eb./N0))/2; %QPSK�����������
semilogy(SNR,BER_true,'-y');
hold on

legend('16QAM simulated','QPSK simulated','QPSK theoretical');
