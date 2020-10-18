clc;clear all;close all;
N = 10000000;
s = source(N); %��Դ���������и���ΪN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%16QAM

Eb = 2.5;%16QAMÿ����������
mu = 0;
SNR = 0 : 5 : 70;
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
    
%     h1 = normrnd(0,sqrt(1/2),N/2,N/4);              %����������������
%     h_i = h1(1:N/4,:);h_q = h1(N/4+1:N/2 ,:);
%     H = h_i + 1i*h_q;
%     s_r = H*(s_c1 + 1i*s_s1).';
% 
%     r1 = s_r + (n_c + 1i*n_s).';
%     
%     W = inv(H'*H)*(H');
%     
%     W1 = W*H;
%     r_ZF = W * r1;             %���������ź�
%     r = zeros(1,N/2);
%     for j= 1:size(r_ZF,1)
%         r(j) = r_ZF(j,j);
%     end
    
    h1 = normrnd(0,sqrt(1/2),40,N/4);              %����������������
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
    
%     W1 = zeros(20,N/4);
%     for j = 0 : N/80-1
%         h = H(:,20*j+1:20*(j+1));
%         W1(:,20*j+1:20*(j+1)) = W(:,20*j+1:20*(j+1))*h;   %���������Ƿ�Ϊ��λ����
%     end
    
    r_ZF = zeros(1,N/4);
    for j = 0 : N/80-1
       r_ZF(20*j+1:20*(j+1)) = W(:,20*j+1:20*(j+1))*(r1(20*j+1:20*(j+1))).';   %���������ź�
    end
    r_c = real(r_ZF);r_s = imag(r_ZF);  
    
    
    
    y = judgement_16QAM(r_c,r_s);     %%16QAM���룬�о����
    BER(i) = error_rate(s,y);        %%���������
end


semilogy(SNR,BER,'b*');
axis([0 60 10^-6 1]);
hold on;
grid on;
xlabel('SNR/dB');ylabel('BER');
title('BER-SNR,Rayleigh');

r = 4*Eb./N0;
BER_true = (3/2*(1 - sqrt(r./(10 + r))) - 9/16*(1 - sqrt(r./(10 + r)).*(4/pi*atan(sqrt((10+r)./r)))))/4; %16QAM�����������
semilogy(SNR,BER_true,'-m');
hold on


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
    
    h1 = normrnd(0,sqrt(1/2),40,N/2);              %����������������
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
    
%     W1 = zeros(20,N/2);
%     for j = 0 : N/40-1
%         h = H(:,20*j+1:20*(j+1));
%         W1(:,20*j+1:20*(j+1)) = W(:,20*j+1:20*(j+1))*h;   %���������Ƿ�Ϊ��λ����
%     end
    
    r_ZF = zeros(1,N/2);
    for j = 0 : N/40-1
       r_ZF(20*j+1:20*(j+1)) = W(:,20*j+1:20*(j+1))*(r1(20*j+1:20*(j+1))).';   %���������ź�
    end
    
    r_c = real(r_ZF);r_s = imag(r_ZF);
    
    y = judgement_QPSK(r_c,r_s);     %%QPSK���룬�о����
    BER(i) = error_rate(s,y);        %%���������
end


semilogy(SNR,BER,'rs');
hold on;

r = 2*Eb./N0;
BER_true = (1 - sqrt(r./(2 + r)) - 1/4*(1 - sqrt(r./(2 + r)).*(4/pi*atan(sqrt((2+r)./r)))))/2;%QPSK�����������
semilogy(SNR,BER_true,'-y');
hold on

legend('16QAM simulated','16QAM theoretical','QPSK simulated','QPSK theoretical');

