%Rayleigh信道 QPSK调制 biterr symerr比较
clear;
nSamp=100;
numSymb=5000;
ts=1/(numSymb*nSamp);
t=(0:numSymb*nSamp-1)*ts;
fd=100;
M=4;
% snr=0:12; 
snr=0:12;
grayencode=[0 1 3 2];
for ii=1:length(snr)
    msg=randsrc(1,numSymb,0:3);%产生发送符号
    msg_gr=grayencode(msg+1);%Gray编码映射
    msg_tx=pskmod(msg_gr,M);%QPSK调制
    msg_tx=rectpulse(msg_tx,nSamp);%矩形脉冲成形
    h=z_rayleigh(8,fd,t)/2;
    msg_tx=h.*msg_tx;
%     h=rayleighchan(ts,10);
%     msg_tx=filter(h,msg_tx);
    msg_rx=awgn(msg_tx,snr(ii));%通过AWGN信道
    msg_rx_down=intdump(msg_rx,nSamp);%匹配滤波相干解调
    msg_gr_demod=pskdemod(msg_rx_down,M); %QPSK解调
    [dummy,graydecod]=sort(grayencode);
    graydecod=graydecod-1;
    msg_demod=graydecod(msg_gr_demod+1);%Gray逆编码
    [errorBit,ber(ii)]=biterr(msg,msg_demod,log2(M));%计算biterr
    [errorSym,ser(ii)]=symerr(msg,msg_demod);%计算symerr
end
% scatterplot(msg_tx(1:100));%发射信号星座图
% scatterplot(msg_rx(1:100));
figure;
r=10.^(snr/10);
w_ser=2*qfunc(sqrt(2*r)*sin(pi/4));
w_ber=1/log2(4)*w_ser;
semilogy(snr,ber,'-ro',snr,ser,'-r*');
% semilogy(snr,ber,'-ro',snr,ser,'-r*',snr,w_ber,'-bo',snr,w_ser,'-b*');
% legend('仿真误比特率','仿真误符号率','理想误比特率','理想误符号率');
% xlabel('信噪比/dB');
% ylabel('ber/ser');
