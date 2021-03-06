%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Indoor
d = 0 : 5 : 70;    %单位km
N_floor = 0 : 10 : 100;    %传播过程中经过楼层的数量
PL = 37 + 30*log10(d) + 18.3*power(N_floor,((N_floor+2)/(N_floor+1)-0.46));

%Pedestrain
%d = 10;
f_c = 2000;  %IMT-2000的载波频率为2000MHz
PL = 40*log10(d) + 30*log10(f_c) + 49;
%Vehicular
%d = 10;
f_c = 2000;
h_b = 30; % 建筑物平均顶层高度向上测量的发射天线高度,单位为m
PL = 40*(1-4*10^-3*h_b)*log10(d) - 18*log10(h_b) + 21*log10(f_c) +80;
