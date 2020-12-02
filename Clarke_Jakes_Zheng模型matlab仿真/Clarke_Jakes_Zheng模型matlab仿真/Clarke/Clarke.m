clear;clc
% 
% 功能：使用成型滤波法仿真Clarke模型
% 作者：hhuyudi
% 参考：MIMO-OFDM Wireless Communications with MATLAB
%        参考书中 FWGN 小节 
%
%最大多普勒频为926Hz，仿真次数为50us，仿真次数为十万次
fd=926;     % 最大Doppler
ts=5e-5;    % 仿真最小时间
fs=1/ts;    % sample frequency
Ns=1e5;     % 仿真次数
M=2^12;     % 做FFT的总点数，可以使用自相关函数通过FFT计算功率谱密度，请自行完成
t=(0:Ns-1)*ts;% 仿真总时间，作图使用
tt=(0:M-1)*ts;% 仿真总时间，作图使用
ff=(-M/2:M/2-1)/(M*ts*fd);% 频率归一化，功率谱密度作图使用
temp=zeros(3,Ns); % 3*Ns的矩阵
                  % temp（1，：）存储仿真结果的自相关函数
                  % temp（2，：）存储ones（1，Ns）的自相关函数，用于归一化
                  % temp（3，：）存储仿真结果其分量的互相关函数
                  
%generat channel information
[h,Nfft,Nifft,Doppler_coeff]=Clarke_model(fd,fs,Ns); % 归一化的Clarke模型，功率为1

%自相关函数，分量的互相关函数，为了结果更准确，这里50次仿真取平均值
for ii=1:50
    [h,Nfft,Nifft,Doppler_coeff]=Clarke_model(fd,fs,Ns);
    yy=xcorr(h);
    yy_cs=xcorr(real(h),imag(h));
    temp(1,:)=temp(1,:)+yy(Ns:length(yy));
    temp(3,:)=temp(3,:)+yy_cs(Ns:length(yy_cs));
end

%画图
%Clarke模型复包络的幅度增益、幅度分布、相位分布
figure(1);
subplot(211)
plot((1:Ns)*ts,10*log10(abs(h)));
axis([0 0.01 -20 10]);
xlabel('时间/s');ylabel('幅度/db');title('信道时域特性');
str=sprintf('channel model by Clarke with fm=%dHz,Ts=%dus',fd,ts*1e6);
title(str);
% 幅度分布
subplot(223)
[f,xi]=ksdensity(abs(h));
plot(xi,f);
cs2=var(h)/2;  % 瑞利分布的方差是实部或虚部的方差，是h方差的一半
r=linspace(0,2,1000);
fx2=r./(cs2).*exp(-r.^2/2/(cs2));
hold on;plot(r,fx2,'r:');hold off;
xlabel('幅度');ylabel('统计次数');title('幅度分布');axis([0 2.5 0 1.0]);
% 相位分布
subplot(224)
[f,xi]=ksdensity(angle(h));
plot(xi,f);
hold on;plot([-pi pi],[1/2/pi 1/2/pi],'r:');hold off;
xlabel('相位/rad');ylabel('统计次数');title('相位分布');axis([-pi pi 0 0.2]);


%%自相关函数和互相关函数计算，功率谱密度请自行编程
temp(1,:)=temp(1,:)/50;
temp(3,:)=temp(3,:)/50;
%用于归一化
yyy=xcorr(ones(1,Ns));
temp(2,:)=yyy(Ns:length(yy));
%自相关函数和分量互相关函数，取M个点，进行归一化
for k=1:M
    simulated_corr(k)=real(temp(1,k))/temp(2,k);
    simulated_corr_cs(k)=real(temp(3,k))/temp(2,k);
end
classical_corr=besselj(0,2*pi*fd*tt);

%%画图：自相关函数和互相关函数
figure(2);
subplot(211)
plot(tt,simulated_corr,'k-');
title('自相关函数'); grid on; 
xlabel('时间差/s');ylabel('相关系数');axis([0 0.004 -0.5 1]);
legend('仿真结果');
subplot(212)
plot(tt,simulated_corr_cs,'k-');
title('互相关函数'); grid on; 
xlabel('时间差/s');ylabel('相关系数');axis([0 0.004 -0.2 0.2]);