% �Ľ�Jakesģ�ͣ�Ҳ����Zhengģ��
clear;clc;
fd=926;
ts=5e-5;
Ns=1e5;
M=2^12;

t=(0:Ns-1)*ts;
tt=(0:M-1)*ts;
ff=[-M/2:M/2-1]/(M*ts*fd);
temp=zeros(3,Ns);

for ii=1:50
    h=z_rayleigh(8,fd,t);
    yy=xcorr(h);
    yy_cs=xcorr(real(h),imag(h));
    temp(1,:)=temp(1,:)+yy(Ns:length(yy));
    temp(3,:)=temp(3,:)+yy_cs(Ns:length(yy_cs));
end

%%����غ����͹������ܶ�
temp(1,:)=temp(1,:)/50;
temp(3,:)=temp(3,:)/50;
%���ڹ�һ��
yyy=xcorr(ones(1,Ns));
temp(2,:)=yyy(Ns:length(yy));
%����غͷ��������
for k=1:M
    simulated_corr(k)=real(temp(1,k))/temp(2,k);
    simulated_corr_cs(k)=real(temp(3,k))/temp(2,k);
end
classical_corr=besselj(0,2*pi*fd*tt);
%�������ܶȣ�������ͼ
classical_y=fftshift(fft(classical_corr));
simulated_y=fftshift(fft(simulated_corr));

%%��ͼ������غ����ͻ���غ���
figure(2);
subplot(211)
plot(tt,classical_corr,'k-',tt,simulated_corr,'r-');
title('����غ���'); grid on; 
xlabel('ʱ���/s');ylabel('���ϵ��');axis([0 0.004 -0.5 1]);
legend('��������','������');
subplot(212)
plot(tt,simulated_corr_cs,'r-');
title('����غ���'); grid on; 
xlabel('ʱ���/s');ylabel('���ϵ��');axis([0 0.004 -0.2 0.2]);



%���Ⱥ���λ�ֲ�����
[f,xi]=ksdensity(abs(h));
[f2,xii]=ksdensity(angle(h));
%rayleigh�ֲ���pdf
cs2=var(h)/2;  %����Ϊʵ�����鲿�ķ���
r=linspace(0,2,1000);
fx2=r./(cs2).*exp(-r.^2/2/(cs2));

%%��ͼ�����ȡ����ȷֲ�����λ�ֲ�
figure(1);
%��������
subplot(211);
plot(t,10*log10(abs(h)));
ylabel('����/dB');xlabel('ʱ��/s');axis([0 0.01 -20 10]);
title('�Ľ�Jakesģ�ͣ�fd=926Hz��ts=50us');
subplot(223);
plot(xi,f);
hold on;plot(r,fx2,'r:');hold off;
ylabel('���ȷֲ�');xlabel('����');title('���ȷֲ�');axis([0 2.5 0 1]);
subplot(224);
plot(xii,f2);
hold on;plot([-pi pi],[1/2/pi 1/2/pi],'r:');hold off;
ylabel('��λ�ֲ�');xlabel('��λ/rad');title('��λ�ֲ�');axis([-pi pi 0 0.2]);