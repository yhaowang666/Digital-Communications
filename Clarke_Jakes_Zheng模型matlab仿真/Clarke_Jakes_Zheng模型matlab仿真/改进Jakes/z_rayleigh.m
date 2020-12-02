function [ h ] = z_rayleigh( M,fd,t )
%Z_RAYLEIGH �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

%MΪ���Ҳ�����fdΪ��������Ƶ�ƣ�tΪʱ������
N=4*M;
wd=2*pi*fd;
zc=zeros(1,length(t));
zs=zeros(1,length(t));
P_nor=sqrt(1/M);
%��������Ƶ�ƺͱ���
for ii=1:M
    phi=2*pi*rand(1,1)-pi;
    psi=2*pi*rand(1,1)-pi;
    theta=2*pi*rand(1,1)-pi;
    alpha(ii)=(2*pi*ii-pi+theta)/4/M;
    zc=zc+cos(cos(alpha(ii))*wd*t+phi);
    zs=zs+cos(sin(alpha(ii))*wd*t+psi);
end
h=P_nor*(zc+1i*zs);
end

