function [h,Nfft,Nifft,doppler_coeff]=Clarke_model(fm,fs,N)
% Clarke model
% function��
%       ʹ�ó����˲�������Clarkeģ�͵ķ���
% inputs:
%       fm   ��������Ƶ��
%       fs   ����Ƶ��
%       N    ������
% output:
%       h    �����ŵ�
% make sure the Nfft=2^n
Nfft=2^max(3,nextpow2(2*fm/fs*N));
Nifft=ceil(Nfft*fs/(2*fm));

% �����������
GI=randn(1,Nfft);
GQ=randn(1,Nfft);
% �����ź�ȡFFT
CGI=fft(GI);
CGQ=fft(GQ);
% ����Nfft�������ղ���ֵ
doppler_coeff=Doppler(fm,Nfft);
% Ƶ����˽����˲�
f_CGI=CGI.*sqrt(doppler_coeff); 
f_CGQ=CGQ.*sqrt(doppler_coeff);

%Ϊ��ʵ��IFFT��ͨ����Nifft-Nfft�����������������������ߴ�
tzeros=zeros(1,Nifft-Nfft);
filtered_CGI=[f_CGI(1:Nfft/2) tzeros f_CGI(Nfft/2+1:Nfft)];
filtered_CGQ=[f_CGQ(1:Nfft/2) tzeros f_CGQ(Nfft/2+1:Nfft)];
hI=ifft(filtered_CGI);
hQ=ifft(filtered_CGQ);

% Clarke�����磺ʵ���鲿�Ƕ���ͬ�ֲ�ͬ����ĸ�˹�ֲ���Ϊ����Ӷ�����Ƶ�ƣ��������˲�����
rayEnvelope=sqrt(abs(hI).^2+abs(hQ).^2);
rayRMS=sqrt(mean(rayEnvelope(1:N).*rayEnvelope(1:N))); % ����Ϊ�˹�һ��
h=complex(real(hI(1:N)),-real(hQ(1:N)))/rayRMS; % h ��Ϊ�ŵ�������
