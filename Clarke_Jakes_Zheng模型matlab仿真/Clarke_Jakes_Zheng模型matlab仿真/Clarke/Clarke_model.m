function [h,Nfft,Nifft,doppler_coeff]=Clarke_model(fm,fs,N)
% Clarke model
% function：
%       使用成型滤波法进行Clarke模型的仿真
% inputs:
%       fm   最大多普勒频移
%       fs   采样频率
%       N    采样数
% output:
%       h    复合信道
% make sure the Nfft=2^n
Nfft=2^max(3,nextpow2(2*fm/fs*N));
Nifft=ceil(Nfft*fs/(2*fm));

% 产生随机过程
GI=randn(1,Nfft);
GQ=randn(1,Nfft);
% 对是信号取FFT
CGI=fft(GI);
CGQ=fft(GQ);
% 生成Nfft个多普勒采样值
doppler_coeff=Doppler(fm,Nfft);
% 频域相乘进行滤波
f_CGI=CGI.*sqrt(doppler_coeff); 
f_CGQ=CGQ.*sqrt(doppler_coeff);

%为了实现IFFT，通过（Nifft-Nfft）个采样补零来调整采样尺寸
tzeros=zeros(1,Nifft-Nfft);
filtered_CGI=[f_CGI(1:Nfft/2) tzeros f_CGI(Nfft/2+1:Nfft)];
filtered_CGQ=[f_CGQ(1:Nfft/2) tzeros f_CGQ(Nfft/2+1:Nfft)];
hI=ifft(filtered_CGI);
hQ=ifft(filtered_CGQ);

% Clarke复包络：实、虚部是独立同分布同方差的高斯分布，为了添加多普勒频移，进行了滤波处理
rayEnvelope=sqrt(abs(hI).^2+abs(hQ).^2);
rayRMS=sqrt(mean(rayEnvelope(1:N).*rayEnvelope(1:N))); % 这是为了归一化
h=complex(real(hI(1:N)),-real(hQ(1:N)))/rayRMS; % h 即为信道复包络
