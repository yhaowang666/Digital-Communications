clc;clear all;close all;
a=0.01:0.01:0.5;N0=2*power(a,2);
N=[2000 20000 200000];
for i=1:3
    s=source(N(i));
    E=0;
    for m=1:N(i)
        E=E+s(m);
    end
    Eb=E/N(i);            %%��ƽ��ÿ��������
    for b=1:length(a)
        SNR(b)=10*log(Eb/N0(b));   %%�������
        n=normrnd(0,a(b),[2,N(i)/2]);   %����˫·����
        n_c=n(1,:);n_s=n(2,:);
        s1_c=zeros(1,N(i)/2);s1_s=zeros(1,N(i)/2);
        for c=1:N(i)/2
            s1_c(c)=s(2*c-1);
            s1_s(c)=s(2*c);
        end                     %%����Դ�ֽ��˫·�ź�
        [s_c,s_s]=QPSK(s1_c,s1_s);     %%����QPSK����
        r_c=s_c+n_c;r_s=s_s+n_s;
        y=judgement_QPSK2(r_c,r_s);     %%�о����
        BER(b)=error_rate(s,y); %%���������
        BER_true(b)=erfc(sqrt(Eb/N0(b)))/2;
    end
    figure(i)
    semilogy(SNR,BER_true,'r');xlabel('SNR/dB');ylabel('SER');
    hold on
    semilogy(SNR,BER,'b')
    hold on
end
