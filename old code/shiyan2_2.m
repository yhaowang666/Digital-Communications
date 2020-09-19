clc;clear all;close all;
a=0.2:0.01:0.4;N0=2*power(a,2);
N=[3000 30000 300000];
for i=1:3
    s=source(N(i));
    E=0;
    Eb=1/3;            %%ƽ��ÿ��������
    for b=1:length(a)
        SNR(b)=10*log(Eb/N0(b));   %%�������
        n=normrnd(0,a(b),[2,N(i)/3]);   %����˫·����
        n_c=n(1,:);n_s=n(2,:);
        s1=zeros(1,N(i)/3);s2=zeros(1,N(i)/3);s3=s1;
        for c=1:N(i)/3
            s1(c)=s(3*c-2);
            s2(c)=s(3*c-1);
            s3(c)=s(3*c);
        end                     %%����Դ�ֽ����·�ź�
        [s_c,s_s]=Eight_PSK(s1,s2,s3);     %%����8PSK����
        r_c=s_c+n_c;r_s=s_s+n_s;
        y=judgement_EightPSK(r_c,r_s);     %%�о����
        SER(b)=S_ER(s,y); %%���������
        SER_true(b)=erfc(sqrt(3*Eb/N0(b))*sin(pi/8));  %%�����������
    end
    figure(i)
    semilogy(SNR,SER_true,'r');xlabel('SNR/dB');ylabel('SER');
    hold on
    semilogy(SNR,SER,'b')
    hold on
end
