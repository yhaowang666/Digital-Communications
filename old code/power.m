SNR=-5:0.1:20;          
N=600000;
s=source(N);   %��Դ����

%QPSK����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eb=1/2;    %%ƽ��ÿ��������
N0=Eb./(power(10,SNR/10));
a=sqrt(N0/2);
    for b=1:length(a)
        n=normrnd(0,a(b),[2,N/2]);   %����˫·����
        n_c=n(1,:);n_s=n(2,:);
        s1_c=zeros(1,N/2);s1_s=zeros(1,N/2);
        for c=1:N/2
            s1_c(c)=s(2*c-1);
            s1_s(c)=s(2*c);
        end                     %%����Դ�ֽ��˫·�ź�
        [s_c,s_s]=QPSK(s1_c,s1_s);     %%����QPSK����
        r_c=s_c+n_c;r_s=s_s+n_s;
        y=judgement_QPSK2(r_c,r_s);     %%�о����
        ber(b)=error_rate(s,y); %%���������
        %BER_true(b)=erfc(sqrt(Eb/N0(b)))/2;
    end
    %semilogy(SNR,BER_true*2,'r');
    semilogy(SNR,ber,'b')
    xlabel('SNR/dB');ylabel('BER');
    hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%8PSK����   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eb=1/3;    %%ƽ��ÿ��������
N0=Eb./(power(10,SNR/10));
a=sqrt(N0/2);
    for b=1:length(a)
        n=normrnd(0,a(b),[2,N/3]);   %����˫·����
        n_c=n(1,:);n_s=n(2,:);
        s1=zeros(1,N/3);s2=zeros(1,N/3);s3=s1;
        for c=1:N/3
            s1(c)=s(3*c-2);
            s2(c)=s(3*c-1);
            s3(c)=s(3*c);
        end                     %%����Դ�ֽ����·�ź�
        [s_c,s_s]=Eight_PSK(s1,s2,s3);     %%����8PSK����
        r_c=s_c+n_c;r_s=s_s+n_s;
        y=judgement_EightPSK(r_c,r_s);     %%�о����
        ber(b)=BER(s,y); %%���������
        %SER_true(b)=erfc(sqrt(3*Eb/N0(b))*sin(pi/8));  %%�����������
    end
    %semilogy(SNR,SER_true,'r');xlabel('SNR/dB');ylabel('SER');
    semilogy(SNR,ber,'r')
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%SM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eb=1/4;    %%ƽ��ÿ��������
N0=Eb./(power(10,SNR/10));
a=sqrt(N0/2);
    for b=1:length(a)
        n=normrnd(0,a(b),[2,N/4]);   %����˫·����
        n_c=n(1,:);n_s=n(2,:);
        sm1=zeros(1,N/4);sm2=zeros(1,N/4);sm3=sm1;sm4=sm1;
        for c=1:N/4
            sm1(c)=s(4*c-3);
            sm2(c)=s(4*c-2);
            sm3(c)=s(4*c-1);
            sm4(c)=s(4*c);
        end                     %%����Դ�ֽ����·�ź�
        l=Antenna_mapping(sm1,sm2);
        H1=[0.5,1,1.5,2];       %4��1�ŵ������ŵ����Խϲ�����
        [s_c,s_s]=QPSK(sm3,sm4);     %%����QPSK����
        for i=1:length(l)
         s_c(i)=s_c(i)*H1(l(i));
         s_s(i)=s_s(i)*H1(l(i));
        end
        r_c=s_c+n_c;r_s=s_s+n_s;
        y=ML(r_c,r_s,H1);     %%�����Ȼ�о�����о����
        ber(b)=BER(s,y); %%���������
        %SER_true(b)=erfc(sqrt(3*Eb/N0(b))*sin(pi/8));  %%�����������
    end
    %semilogy(SNR,SER_true,'r');xlabel('SNR/dB');ylabel('SER');
    semilogy(SNR,ber,'k')
    hold on
    
    
Eb=1/4;    %%ƽ��ÿ��������
N0=Eb./(power(10,SNR/10));
a=sqrt(N0/2);
    for b=1:length(a)
        n=normrnd(0,a(b),[2,N/4]);   %����˫·����
        n_c=n(1,:);n_s=n(2,:);
        sm1=zeros(1,N/4);sm2=zeros(1,N/4);sm3=sm1;sm4=sm1;
        for c=1:N/4
            sm1(c)=s(4*c-3);
            sm2(c)=s(4*c-2);
            sm3(c)=s(4*c-1);
            sm4(c)=s(4*c);
        end                     %%����Դ�ֽ����·�ź�
        l=Antenna_mapping(sm1,sm2);
        H2=[2,4,6,8];       %4��1�ŵ������ŵ���������
        [s_c,s_s]=QPSK(sm3,sm4);     %%����QPSK����
        for i=1:length(l)
         s_c(i)=s_c(i)*H2(l(i));
         s_s(i)=s_s(i)*H2(l(i));
        end
        r_c=s_c+n_c;r_s=s_s+n_s;
        y=ML(r_c,r_s,H2);     %%�����Ȼ�о�����о����
        ber(b)=BER(s,y); %%���������
        %SER_true(b)=erfc(sqrt(3*Eb/N0(b))*sin(pi/8));  %%�����������
    end
    %semilogy(SNR,SER_true,'r');xlabel('SNR/dB');ylabel('SER');
    semilogy(SNR,ber,'g')
    hold on
    legend('QPSK','8PSK','4��1SM�������ŵ�','4��1SM�����ŵ�')