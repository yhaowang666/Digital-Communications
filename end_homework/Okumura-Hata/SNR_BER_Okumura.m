clc;clear all;close all;
N = 10000000;
s = source(N); %��Դ���������и���ΪN
SNR = 10 : 10 : 90;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_c = 1500;%��λΪMHZ
d = 4;%��λΪkm 
h_b = 200;%��λΪm 
h_m = 3; %��λΪm 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Okumura-Hataģ�ͣ��ز���ΧΪ150MHZ��1500MHZ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��С�ͳ���
% a = (1.1*log(f_c) - 0.7)*h_m - (1.56*log(f_c) - 0.8);
% C = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %������
% for i = 1:length(f_c)
%     if(f_c(i) < 300)
%         a(i) = 8.29*(log(1.54*h_m))^2 - 1.1;
%     else
%         a(i) = 3.2*(log(11.75*h_m))^2 -4.97;
%     end
% end
% C = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %��������
% a = (1.1*log(f_c) - 0.7)*h_m - (1.56*log(f_c) - 0.8);
% C = -2*(log(f_c/28)).^2 - 5.4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %ũ��
a = (1.1*log(f_c) - 0.7)*h_m - (1.56*log(f_c) - 0.8);
C = -4.78*(log(f_c)).^2 + 18.33*log(f_c) - 40.98;
A = 69.55 + 26.16*log(f_c) - 13.82*log(h_b) - a;
B = 44.9 - 6.55*log(h_b);
PL = A + B*log(d) + C;

EL = power(10,PL/10);    %������ģ����������������������Ϊ��ǰ��1/EL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%16QAM
%��·��˥��
Eb = 2.5;%16QAMÿ����������
mu = 0;
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
    
    [s_c1,s_s1] = QAM(s1);     %����16QAM����
    
    h1 = normrnd(0,sqrt(1/(2*EL)),40,N/4);              %����������������
    h_i = h1(1:20,:);h_q = h1(21:40 ,:);
    H = h_i + 1i*h_q;
    S = (s_c1 + 1i*s_s1).';
    
    for j = 0 : N/80-1                %����ת����ÿһ��20bit����С���㸴�Ӷ�
        s_r(20*j+1:20*(j+1)) = H(:,20*j+1:20*(j+1))*S(20*j+1:20*(j+1));
    end

    r1 = s_r + n_c + 1i*n_s;
    
    W = zeros(20,N/4);
    for j = 0 : N/80-1
        h = H(:,20*j+1:20*(j+1));
        W(:,20*j+1:20*(j+1)) = (h'*h)\h';
    end
    
    r_ZF = zeros(1,N/4);
    for j = 0 : N/80-1
       r_ZF(20*j+1:20*(j+1)) = W(:,20*j+1:20*(j+1))*(r1(20*j+1:20*(j+1))).';   %���������ź�
    end
    r_c = real(r_ZF);r_s = imag(r_ZF);  
    
    y = judgement_16QAM(r_c,r_s);     %%16QAM���룬�о����
    BER(i) = error_rate(s,y);         %%���������
end

semilogy(SNR,BER,'-md');
%axis([0 60 10^-6 1]);
hold on;
grid on;
xlabel('SNR/dB');ylabel('BER');
title({'BER-SNR,Okumura-Hata';['rural area,distance=',num2str(d),'f_c=',num2str(f_c)]});

% %��·��˥��
% for i = 1:length(sigma)
%     n = normrnd(mu,sigma(i),[2,N/4]);   %�������Ӹ�˹�ֲ���˫·����
%     n_c = n(1,:);n_s = n(2,:);
%     s1 = zeros(4,N/4);
% 
%     for c = 1:N/4
%         s1(1,c) = s(4*c-3);
%         s1(2,c) = s(4*c-2);
%         s1(3,c) = s(4*c-1);
%         s1(4,c) = s(4*c);
%     end                     %����Դ�ֽ����·�ź�
%     
%     [s_c1,s_s1] = QAM(s1);     %����16QAM����
%     
%     h1 = normrnd(0,sqrt(1/2),40,N/4);              %����������������
%     h_i = h1(1:20,:);h_q = h1(21:40 ,:);
%     H = h_i + 1i*h_q;
%     S = (s_c1 + 1i*s_s1).';
%     
%     for j = 0 : N/80-1                %����ת����ÿһ��20bit����С���㸴�Ӷ�
%         s_r(20*j+1:20*(j+1)) = H(:,20*j+1:20*(j+1))*S(20*j+1:20*(j+1));
%     end
% 
%     r1 = s_r + n_c + 1i*n_s;
%     
%     W = zeros(20,N/4);
%     for j = 0 : N/80-1
%         h = H(:,20*j+1:20*(j+1));
%         W(:,20*j+1:20*(j+1)) = (h'*h)\h';
%     end
%     
%     r_ZF = zeros(1,N/4);
%     for j = 0 : N/80-1
%        r_ZF(20*j+1:20*(j+1)) = W(:,20*j+1:20*(j+1))*(r1(20*j+1:20*(j+1))).';   %���������ź�
%     end
%     r_c = real(r_ZF);r_s = imag(r_ZF);  
%     
%     y = judgement_16QAM(r_c,r_s);     %%16QAM���룬�о����
%     BER(i) = error_rate(s,y);         %%���������
% end
% semilogy(SNR,BER,'--md');
% hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QPSK
%��·��˥��
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
    
    h1 = normrnd(0,sqrt(1/(2*EL)),40,N/2);              %����������������
    h_i = h1(1:20,:);h_q = h1(21:40 ,:);
    H = h_i + 1i*h_q;
    S = (s_c + 1i*s_s).';
    
    for j = 0 : N/40-1                %����ת����ÿһ��20bit����С���㸴�Ӷ�
        s_r(20*j+1:20*(j+1)) = H(:,20*j+1:20*(j+1))*S(20*j+1:20*(j+1));
    end

    r1 = s_r + n_c + 1i*n_s;
    
    W = zeros(20,N/2);
    for j = 0 : N/40-1
        h = H(:,20*j+1:20*(j+1));
        W(:,20*j+1:20*(j+1)) = (h'*h)\h';
    end
    
    r_ZF = zeros(1,N/2);
    for j = 0 : N/40-1
       r_ZF(20*j+1:20*(j+1)) = W(:,20*j+1:20*(j+1))*(r1(20*j+1:20*(j+1))).';   %���������ź�
    end
    
    r_c = real(r_ZF);r_s = imag(r_ZF);
    
    y = judgement_QPSK(r_c,r_s);     %%QPSK���룬�о����
    BER(i) = error_rate(s,y);        %%���������
end
semilogy(SNR,BER,'-gs');
hold on;

% %��·��˥��
% for i = 1:length(sigma)
%     n = normrnd(mu,sigma(i),[2,N/2]);   %�������Ӹ�˹�ֲ���˫·����
%     n_c = n(1,:);n_s = n(2,:);
%     s1_c = zeros(1,N/2);s1_s = zeros(1,N/2);
% 
%     for c = 1:N/2
%         s1_c(c) = s(2*c-1);
%         s1_s(c) = s(2*c);
%     end                     %����Դ�ֽ��˫·�ź�
%     
%     [s_c,s_s] = QPSK(s1_c,s1_s);     %����QPSK����
%     
%     h1 = normrnd(0,sqrt(1/2),40,N/2);              %����������������
%     h_i = h1(1:20,:);h_q = h1(21:40 ,:);
%     H = h_i + 1i*h_q;
%     S = (s_c + 1i*s_s).';
%     
%     for j = 0 : N/40-1                %����ת����ÿһ��20bit����С���㸴�Ӷ�
%         s_r(20*j+1:20*(j+1)) = H(:,20*j+1:20*(j+1))*S(20*j+1:20*(j+1));
%     end
% 
%     r1 = s_r + n_c + 1i*n_s;
%     
%     W = zeros(20,N/2);
%     for j = 0 : N/40-1
%         h = H(:,20*j+1:20*(j+1));
%         W(:,20*j+1:20*(j+1)) = (h'*h)\h';
%     end
%     
%     r_ZF = zeros(1,N/2);
%     for j = 0 : N/40-1
%        r_ZF(20*j+1:20*(j+1)) = W(:,20*j+1:20*(j+1))*(r1(20*j+1:20*(j+1))).';   %���������ź�
%     end
%     
%     r_c = real(r_ZF);r_s = imag(r_ZF);
%     
%     y = judgement_QPSK(r_c,r_s);     %%QPSK���룬�о����
%     BER(i) = error_rate(s,y);        %%���������
% end
% semilogy(SNR,BER,'--gs');
% hold on;

legend('16QAM-PL','QPSK-PL');
