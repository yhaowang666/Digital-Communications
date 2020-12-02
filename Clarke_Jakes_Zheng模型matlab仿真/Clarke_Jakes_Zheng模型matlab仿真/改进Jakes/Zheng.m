% 改进Jakes模型，也就是Zheng模型
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

%%自相关函数和功率谱密度
temp(1,:)=temp(1,:)/50;
temp(3,:)=temp(3,:)/50;
%用于归一化
yyy=xcorr(ones(1,Ns));
temp(2,:)=yyy(Ns:length(yy));
%自相关和分量互相关
for k=1:M
    simulated_corr(k)=real(temp(1,k))/temp(2,k);
    simulated_corr_cs(k)=real(temp(3,k))/temp(2,k);
end
classical_corr=besselj(0,2*pi*fd*tt);
%功率谱密度，不再作图
classical_y=fftshift(fft(classical_corr));
simulated_y=fftshift(fft(simulated_corr));

%%画图：自相关函数和互相关函数
figure(2);
subplot(211)
plot(tt,classical_corr,'k-',tt,simulated_corr,'r-');
title('自相关函数'); grid on; 
xlabel('时间差/s');ylabel('相关系数');axis([0 0.004 -0.5 1]);
legend('理想特性','仿真结果');
subplot(212)
plot(tt,simulated_corr_cs,'r-');
title('互相关函数'); grid on; 
xlabel('时间差/s');ylabel('相关系数');axis([0 0.004 -0.2 0.2]);



%幅度和相位分布函数
[f,xi]=ksdensity(abs(h));
[f2,xii]=ksdensity(angle(h));
%rayleigh分布的pdf
cs2=var(h)/2;  %方差为实部或虚部的方差
r=linspace(0,2,1000);
fx2=r./(cs2).*exp(-r.^2/2/(cs2));

%%画图：幅度、幅度分布、相位分布
figure(1);
%幅度增益
subplot(211);
plot(t,10*log10(abs(h)));
ylabel('幅度/dB');xlabel('时间/s');axis([0 0.01 -20 10]);
title('改进Jakes模型，fd=926Hz，ts=50us');
subplot(223);
plot(xi,f);
hold on;plot(r,fx2,'r:');hold off;
ylabel('幅度分布');xlabel('幅度');title('幅度分布');axis([0 2.5 0 1]);
subplot(224);
plot(xii,f2);
hold on;plot([-pi pi],[1/2/pi 1/2/pi],'r:');hold off;
ylabel('相位分布');xlabel('相位/rad');title('相位分布');axis([-pi pi 0 0.2]);