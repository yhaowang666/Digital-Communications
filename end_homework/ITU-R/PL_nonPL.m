clc;clear all;close all;
N = 10000000;
s = source(N); %信源产生，序列个数为N
SNR = 20 : 5 : 90;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Indoor
% d = 1;    %单位km
% N_floor = 2;    %传播过程中经过楼层的数量
% PL = 37 + 30*log10(d) + 18.3*power(N_floor,((N_floor+2)/(N_floor+1)-0.46));
% %Pedestrain
% d = 10;
% f_c = 2000;  %IMT-2000的载波频率为2000MHz
% PL = 40*log10(d) + 30*log10(f_c) + 49;
% %Vehicular
% d = 10;
% f_c = 2000;
% h_b = 30; % 建筑物平均顶层高度向上测量的发射天线高度,单位为m
% PL = 40*(1-4*10^-3*h_b)*log10(d) - 18*log10(h_b) + 21*log10(f_c) +80;
PL = 4;
EL = power(10,PL/10);    %能量损耗，产生瑞利乘性噪声方差变为以前的1/EL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%16QAM
%有路径衰落
Eb = 2.5;%16QAM每个比特能量
mu = 0;
BER = zeros(1,length(SNR));
N0 = Eb./(power(10,SNR/10));
sigma = sqrt(N0/2); %计算噪声的标准差

for i = 1:length(sigma)
    n = normrnd(mu,sigma(i),[2,N/4]);   %产生服从高斯分布的双路噪声
    n_c = n(1,:);n_s = n(2,:);
    s1 = zeros(4,N/4);

    for c = 1:N/4
        s1(1,c) = s(4*c-3);
        s1(2,c) = s(4*c-2);
        s1(3,c) = s(4*c-1);
        s1(4,c) = s(4*c);
    end                     %将信源分解成四路信号
    
    [s_c1,s_s1] = QAM(s1);     %进行16QAM编码
    
    h1 = normrnd(0,sqrt(1/(2*EL)),40,N/4);              %产生瑞利乘性噪声
    h_i = h1(1:20,:);h_q = h1(21:40 ,:);
    H = h_i + 1i*h_q;
    S = (s_c1 + 1i*s_s1).';
    
    for j = 0 : N/80-1                %串并转换，每一列20bit，减小计算复杂度
        s_r(20*j+1:20*(j+1)) = H(:,20*j+1:20*(j+1))*S(20*j+1:20*(j+1));
    end

    r1 = s_r + n_c + 1i*n_s;
    
    W = zeros(20,N/4);
    for j = 0 : N/80-1
        h = H(:,20*j+1:20*(j+1));
        W(:,20*j+1:20*(j+1)) = (h'*h)\h';
    end
    
    r_ZF = zeros(1,N/4);
    for j = 0 : N/80-1
       r_ZF(20*j+1:20*(j+1)) = W(:,20*j+1:20*(j+1))*(r1(20*j+1:20*(j+1))).';   %均衡后输出信号
    end
    r_c = real(r_ZF);r_s = imag(r_ZF);  
    
    y = judgement_16QAM(r_c,r_s);     %%16QAM解码，判决输出
    BER(i) = error_rate(s,y);         %%求误比特率
end

semilogy(SNR,BER,'-md');
%axis([0 60 10^-6 1]);
hold on;
grid on;
xlabel('SNR/dB');ylabel('BER');
title(['BER-SNR,ITU-R models,PL = ',num2str(PL)]);

%无路径衰落
for i = 1:length(sigma)
    n = normrnd(mu,sigma(i),[2,N/4]);   %产生服从高斯分布的双路噪声
    n_c = n(1,:);n_s = n(2,:);
    s1 = zeros(4,N/4);

    for c = 1:N/4
        s1(1,c) = s(4*c-3);
        s1(2,c) = s(4*c-2);
        s1(3,c) = s(4*c-1);
        s1(4,c) = s(4*c);
    end                     %将信源分解成四路信号
    
    [s_c1,s_s1] = QAM(s1);     %进行16QAM编码
    
    h1 = normrnd(0,sqrt(1/2),40,N/4);              %产生瑞利乘性噪声
    h_i = h1(1:20,:);h_q = h1(21:40 ,:);
    H = h_i + 1i*h_q;
    S = (s_c1 + 1i*s_s1).';
    
    for j = 0 : N/80-1                %串并转换，每一列20bit，减小计算复杂度
        s_r(20*j+1:20*(j+1)) = H(:,20*j+1:20*(j+1))*S(20*j+1:20*(j+1));
    end

    r1 = s_r + n_c + 1i*n_s;
    
    W = zeros(20,N/4);
    for j = 0 : N/80-1
        h = H(:,20*j+1:20*(j+1));
        W(:,20*j+1:20*(j+1)) = (h'*h)\h';
    end
    
    r_ZF = zeros(1,N/4);
    for j = 0 : N/80-1
       r_ZF(20*j+1:20*(j+1)) = W(:,20*j+1:20*(j+1))*(r1(20*j+1:20*(j+1))).';   %均衡后输出信号
    end
    r_c = real(r_ZF);r_s = imag(r_ZF);  
    
    y = judgement_16QAM(r_c,r_s);     %%16QAM解码，判决输出
    BER(i) = error_rate(s,y);         %%求误比特率
end
semilogy(SNR,BER,'--md');
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QPSK
%有路径衰落
Eb = 1/2;%QPSK每个比特能量为1/2
N0 = Eb./(power(10,SNR/10));
sigma = sqrt(N0/2); %计算噪声的标准差

for i = 1:length(sigma)
    n = normrnd(mu,sigma(i),[2,N/2]);   %产生服从高斯分布的双路噪声
    n_c = n(1,:);n_s = n(2,:);
    s1_c = zeros(1,N/2);s1_s = zeros(1,N/2);

    for c = 1:N/2
        s1_c(c) = s(2*c-1);
        s1_s(c) = s(2*c);
    end                     %将信源分解成双路信号
    
    [s_c,s_s] = QPSK(s1_c,s1_s);     %进行QPSK编码
    
    h1 = normrnd(0,sqrt(1/(2*EL)),40,N/2);              %产生瑞利乘性噪声
    h_i = h1(1:20,:);h_q = h1(21:40 ,:);
    H = h_i + 1i*h_q;
    S = (s_c + 1i*s_s).';
    
    for j = 0 : N/40-1                %串并转换，每一列20bit，减小计算复杂度
        s_r(20*j+1:20*(j+1)) = H(:,20*j+1:20*(j+1))*S(20*j+1:20*(j+1));
    end

    r1 = s_r + n_c + 1i*n_s;
    
    W = zeros(20,N/2);
    for j = 0 : N/40-1
        h = H(:,20*j+1:20*(j+1));
        W(:,20*j+1:20*(j+1)) = (h'*h)\h';
    end
    
    r_ZF = zeros(1,N/2);
    for j = 0 : N/40-1
       r_ZF(20*j+1:20*(j+1)) = W(:,20*j+1:20*(j+1))*(r1(20*j+1:20*(j+1))).';   %均衡后输出信号
    end
    
    r_c = real(r_ZF);r_s = imag(r_ZF);
    
    y = judgement_QPSK(r_c,r_s);     %%QPSK解码，判决输出
    BER(i) = error_rate(s,y);        %%求误比特率
end
semilogy(SNR,BER,'-gs');
hold on;

%无路径衰落
for i = 1:length(sigma)
    n = normrnd(mu,sigma(i),[2,N/2]);   %产生服从高斯分布的双路噪声
    n_c = n(1,:);n_s = n(2,:);
    s1_c = zeros(1,N/2);s1_s = zeros(1,N/2);

    for c = 1:N/2
        s1_c(c) = s(2*c-1);
        s1_s(c) = s(2*c);
    end                     %将信源分解成双路信号
    
    [s_c,s_s] = QPSK(s1_c,s1_s);     %进行QPSK编码
    
    h1 = normrnd(0,sqrt(1/2),40,N/2);              %产生瑞利乘性噪声
    h_i = h1(1:20,:);h_q = h1(21:40 ,:);
    H = h_i + 1i*h_q;
    S = (s_c + 1i*s_s).';
    
    for j = 0 : N/40-1                %串并转换，每一列20bit，减小计算复杂度
        s_r(20*j+1:20*(j+1)) = H(:,20*j+1:20*(j+1))*S(20*j+1:20*(j+1));
    end

    r1 = s_r + n_c + 1i*n_s;
    
    W = zeros(20,N/2);
    for j = 0 : N/40-1
        h = H(:,20*j+1:20*(j+1));
        W(:,20*j+1:20*(j+1)) = (h'*h)\h';
    end
    
    r_ZF = zeros(1,N/2);
    for j = 0 : N/40-1
       r_ZF(20*j+1:20*(j+1)) = W(:,20*j+1:20*(j+1))*(r1(20*j+1:20*(j+1))).';   %均衡后输出信号
    end
    
    r_c = real(r_ZF);r_s = imag(r_ZF);
    
    y = judgement_QPSK(r_c,r_s);     %%QPSK解码，判决输出
    BER(i) = error_rate(s,y);        %%求误比特率
end
semilogy(SNR,BER,'--gs');
hold on;

legend('16QAM-PL','16QAM-nonPL','QPSK-PL','QPSK-nonPL');

