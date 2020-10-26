clc;clear all;close all;
N = 10000000;
s = source(N); %��Դ���������и���ΪN
SNR = 0 : 1 : 12;

s1 = hamming_encoding(s);  %(7,4)��������
N1 = length(s1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QPSK
Eb = 1/2;%QPSKÿ����������Ϊ1/2
N0 = Eb./(power(10,SNR/10));
sigma = sqrt(N0/2); %���������ı�׼��
mu = 0;

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


semilogy(SNR,BER,'-bd');
axis([0 12 10^-6 1])
hold on;
grid on;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QPSK Ӳ�о�
Eb = 7/8;%�����QPSKÿ����������Ϊ1/2�������ǰΪ7/8
N0 = Eb./(power(10,SNR/10));
sigma = sqrt(N0/2); %���������ı�׼��
for i = 1:length(sigma)
    n = normrnd(mu,sigma(i),[2,N1/2]);   %�������Ӹ�˹�ֲ���˫·����
    n_c = n(1,:);n_s = n(2,:);
    s1_c = zeros(1,N1/2);s1_s = zeros(1,N1/2);

    for c = 1:N1/2
        s1_c(c) = s1(2*c-1);
        s1_s(c) = s1(2*c);
    end                     %����Դ�ֽ��˫·�ź�
    
    [s_c,s_s] = QPSK(s1_c,s1_s);     %����QPSK����
    r_c = s_c + n_c;r_s = s_s + n_s;
    
    y1 = judgement_QPSK(r_c,r_s);     %%QPSK���룬�о������Ӳ�о���ֱ�����0��1
    
    y = hamming_decoding(y1);         %��С������������
    
    BER(i) = error_rate(s,y);         %���������
end

semilogy(SNR,BER,'-rs');
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QPSK ���о�
% for i = 1:length(sigma)
%     n = normrnd(mu,sigma(i),[2,N1/2]);   %�������Ӹ�˹�ֲ���˫·����
%     n_c = n(1,:);n_s = n(2,:);
%     s1_c = zeros(1,N1/2);s1_s = zeros(1,N1/2);
% 
%     for c = 1:N1/2
%         s1_c(c) = s1(2*c-1);
%         s1_s(c) = s1(2*c);
%     end                     %����Դ�ֽ��˫·�ź�
%     
%     [s_c,s_s] = QPSK(s1_c,s1_s);     %����QPSK����
%     r_c = s_c + n_c;r_s = s_s + n_s;
%     
%     y1 = judgement_QPSK(r_c,r_s);     %%QPSK���룬�о������Ӳ�о���ֱ�����0��1
%     
%     y = hamming_decoding2(y1);         %��С������������
%     
%     BER(i) = error_rate(s,y);        %%���������
% end
for i = 1:length(sigma)
    n = normrnd(mu,sigma(i),[2,N1/2]);   %�������Ӹ�˹�ֲ���˫·����
    n_c = n(1,:);n_s = n(2,:);
    s1_c = zeros(1,N1/2);s1_s = zeros(1,N1/2);

    for c = 1:N1/2
        s1_c(c) = s1(2*c-1);
        s1_s(c) = s1(2*c);
    end                     %����Դ�ֽ��˫·�ź�
    
    [s_c,s_s] = QPSK(s1_c,s1_s);     %����QPSK����
    r_c = s_c + n_c;r_s = s_s + n_s;
    
    for j = 1 : length(r_c)
        r(j*2-1) = r_c(j);
        r(j*2) = r_s(j);
    end
    y = hamming_decoding1(r);         %��С������������
    
    BER(i) = error_rate(s,y);        %%���������
end

semilogy(SNR,BER,'-kp');
hold on;
















legend('Uncoded','Hard decision','Soft decision');

