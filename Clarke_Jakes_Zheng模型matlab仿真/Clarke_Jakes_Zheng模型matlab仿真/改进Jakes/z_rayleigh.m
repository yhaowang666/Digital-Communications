function [ h ] = z_rayleigh( M,fd,t )
%Z_RAYLEIGH 此处显示有关此函数的摘要
%   此处显示详细说明

%M为正弦波数，fd为最大多普勒频移，t为时间序列
N=4*M;
wd=2*pi*fd;
zc=zeros(1,length(t));
zs=zeros(1,length(t));
P_nor=sqrt(1/M);
%引入的随机频移和变量
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

