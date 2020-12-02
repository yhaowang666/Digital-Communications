% plot Jakes channel model
clear;clc
fd=926;     % 最大Doppler
ts_mu=50;   % 持续时间
scale=1e-6;
ts=ts_mu*scale; %采样时间
fs=1/ts;    % sample frequency
Ns=1e5;     % 采样数

M=2^12;
t=(0:Ns-1)*ts;
tt=(0:M-1)*ts;
ff=[-M/2:M/2-1]/(M*ts*fd);
temp=zeros(3,Ns);

%generat channel information
for ii=1:50
    t_state=0;  % 开始采样时刻
    [h,t_state]=Jakes_model(fd,ts,Ns,t_state,1,0);  % generate channel
    yy=xcorr(h);
    yy_cs=xcorr(real(h),imag(h));
    temp(1,:)=temp(1,:)+yy(Ns:length(yy));
    temp(3,:)=temp(3,:)+yy_cs(Ns:length(yy_cs));
end

figure(1);
subplot(211)
plot((1:Ns)*ts,10*log10(abs(h)));
%axis([0 0.5 -20 5]);
axis([0 0.01 -20 10]);
xlabel('时间/s');ylabel('幅度/db');title('信道时域特性');
str=sprintf('channel model by Jakes with fm=%d[Hz],Ts=%d[us]',fd,ts_mu);
title(str);
% 信道包络
subplot(223)
[f,xi]=ksdensity(abs(h));
plot(xi,f);
cs2=var(h)/2;  %公式乘上了系数sqrt(2)
r=linspace(0,2,1000);
fx2=r./(cs2).*exp(-r.^2/2/(cs2));
hold on;plot(r,fx2,'r:');hold off;
xlabel('幅度');ylabel('统计次数');title('幅度分布');axis([0 2.5 0 1.0]);
% 信道相位
subplot(224)
[f,xi]=ksdensity(angle(h));
plot(xi,f);
hold on;plot([-pi pi],[1/2/pi 1/2/pi],'r:');hold off;
xlabel('相位/rad');ylabel('统计次数');title('相位分布');axis([-pi pi 0 0.2]);


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
xlabel('时间差/s');ylabel('相关系数');axis([0 0.004 -0.5 0.5]);