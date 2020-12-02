clc;clear all;close all
N = 1000;
s = source(N); %��Դ���������и���ΪN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%16QAM

Es = 10;%16QAMÿ����������
mu = 0;
SNR = 0 :2 :30;
SER = zeros(1,length(SNR));
N0 = Es./(power(10,SNR/10));
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
    
    h1 = normrnd(0,sqrt(1/2),N/2,N/4);              %����������������
    h_i = h1(1:N/4,:);h_q = h1(N/4+1:N/2 ,:);
    H = h_i + 1i*h_q;
    s_r = H*(s_c + 1i*s_s).';

    r1 = s_r + n_c + 1i*n_s;
    
    W = inv(H'*H)*(H');
    
    W1 = W*H;
    r_ZF = W * r1;             %���������ź�
    r = zeros(1,N/2);
    for j= 1:size(r_ZF,1)
        r(j) = r_ZF(j,j);
    end
    
    r_c = real(r);r_s = imag(r);   
    
    y = judgement_16QAM(r_c,r_s);     %%16QAM���룬�о����
    SER(i) = symbol_error_16QAM(s,y);        %%���������
end

figure(2)
semilogy(SNR,SER,'b*');
%axis([0 20 10^-6 1]);
hold on;
grid on;
xlabel('SNR/dB');ylabel('BER');
title('BER-SNR,AWGN');

% SER_true = 3*qfunc(sqrt(Es/5./N0)) - 9/4*(qfunc(sqrt(Es/5./N0)).^2); %16QAM�����������
% semilogy(SNR,SER_true,'-m');
% hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QPSK
Es = 1;%QPSKÿ����������Ϊ1
N0 = Es./(power(10,SNR/10));
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

%     h = normrnd(0,sqrt(1/2),2,N/2);              %����������������
%     h_i = h(1,:);h_q = h(2,:);
%     
%     H = h_i + 1i*h_q;     
%     W = inv(H'*H )*(H');   %�������
%     n_ZF = W * (n_c + 1i*n_s);  
%     n = zeros(1,N/2);
%     for j= 1:size(n_ZF,1)
%         n(j) = n_ZF(j,j);
%     end
%     
%     r = s_c + 1i*s_s + n;            %���������ź�
%     
%     r_c = real(r);r_s = imag(r);

    h1 = normrnd(0,sqrt(1/2),N,N/2);              %����������������
    h_i = h1(1:N/2,:);h_q = h1(N/2+1:N ,:);
    H = h_i + 1i*h_q;
    S = (s_c + 1i*s_s).';
    s_r = H*S;
    r1 = s_r + n_c + 1i*n_s;
    
    W = inv(H'*H)*(H');
    
    W1 = W*H;
    r_ZF = W * r1;             %���������ź�
    r = zeros(1,N/2);
    for j= 1:size(r_ZF,1)
        r(j) = r_ZF(j,j);
    end
    
    r_c = real(r);r_s = imag(r);

    
    y = judgement_QPSK(r_c,r_s);     %%QPSK���룬�о����
    SER(i) = symbol_error_QPSK(s,y);        %%���������
end

semilogy(SNR,SER,'rs');
hold on;

% SER_true = erfc(sqrt(Es./N0/2)); %QPSK�����������
% semilogy(SNR,SER_true,'-y');
% hold on
% 
% legend('16QAM simulated','16QAM theoretical','QPSK simulated','QPSK theoretical');