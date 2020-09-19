clc;clear all;close all;

N=1000;
s = source(N); %��Դ���������и���ΪN
Eb = 1;
mu = 0;
SNR = 15;
N0 = Eb./(power(10,SNR/10));
sigma = sqrt(N0/2); %���������ı�׼��



for i =1:length(sigma)
    n = normrnd(sigma(i),2,N/2);   %�������������ֲ���˫·����
    n_c=n(1,:);n_s=n(2,:);
    s1_c=zeros(1,N/2);s1_s=zeros(1,N/2);

    for c=1:N/2
        s1_c(c)=s(2*c-1);
        s1_s(c)=s(2*c);
    end                     %����Դ�ֽ��˫·�ź�
    
    [s_c,s_s] = QPSK(s1_c,s1_s);     %����QPSK����
    r_c = s_c + n_c;r_s = s_s + n_s;
    figure(i)
    scatter(r_c,r_s)
    title(sprintf('SNR = %d',SNR(i)));
end

