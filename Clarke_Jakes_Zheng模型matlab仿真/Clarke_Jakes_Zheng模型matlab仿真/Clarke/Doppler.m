% ��������Ƶ��Ϊfdʱ��Ƶ��
function y=Doppler(fd,Nfft)
% fd ��������Ƶ��
% Nfft Ƶ�����
df=2*fd/Nfft;
%��һ��DCΪ0
f(1)=0;
y(1)=1.5/(pi*fd);
%�����׵���������
for i=2:Nfft/2
    f(i)=(i-1)*df;
    y([i Nfft-i+2])=1.5/(pi*fd*sqrt(1-(f(i)/fd)^2));
end
% ʹ�����3��������ʵ�ֶ���ʽ���
% ���fft��fftshift�������������������ڸ�ɶ
nFitPoints=3;
kk=(Nfft/2-nFitPoints:Nfft/2);
polyFreq=polyfit(f(kk),y(kk),nFitPoints);
y((Nfft/2)+1)=polyval(polyFreq,f(Nfft/2)+df);
