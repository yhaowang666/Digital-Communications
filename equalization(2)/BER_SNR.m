clc;clear all;close all;
N = 10000;
s = source(N); %��Դ���������и���ΪN
h1 = [0.1 -0.3 0.6 0.8 1 -0.9 0.5 0.2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%16QAM

Eb = (2*sqrt(2)+sqrt(10))/4;%16QAMÿ����������
mu = 0;
SNR = -5 :1 : 20;
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


semilogy(SNR,BER,'-b*');
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
    
    [s_i,s_q] = QPSK(s1_c,s1_s);     %����QPSK����
    s_c1 = conv(h1,s_i);
    s_s1 = conv(h1,s_q);             %�����ྶ�ŵ�
    
    h = normrnd(0,sqrt(1/2),2,N/2);              %����������������
    h_i = h(1,:);h_q = h(2,:);
    s_c = s_c1.*h_i - s_s1.*h_q ;s_s = s_c1.*h_q + s_s1.*h_i;
    
    r_c = s_c + n_c;r_s = s_s + n_s;
    
<<<<<<< HEAD
    c = force_zero(h1 ,5);
    r_c1 = conv(c,r_c);
    r_s1 = conv(c,r_s);
    
    y = judgement_QPSK(r_c1,r_s1);     %%QPSK���룬�о����

    h1 = normrnd(0,sqrt(1/2),1,7);
    h1 = h1 ./max(h1);
    n = 5; %��ͷϵ������Ϊ2*n
    c = force_zero(h1,5);
    r_c1 = conv(c,r_c);r_s1 = conv(c,r_s);
    y1 = judgement_QPSK(r_c1,r_s1);
    y = y1(1+2*n:N+2*n);%%QPSK���룬�о����

    BER(i) = error_rate(s,y);        %%���������
end


semilogy(SNR,BER,'-rs');
hold on;
legend('16QAM','QPSK');