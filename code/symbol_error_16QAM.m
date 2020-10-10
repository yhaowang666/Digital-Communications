function SER = symbol_error_16QAM(x,y)
N = length(x)/4;
ser = 0;
for i = 1:N
    if x(4*i-3)~= y(4*i-3) || x(4*i-2)~= y(4*i-2) || x(4*i-1)~= y(4*i-1) || x(4*i)~= y(4*i)
        ser =ser+1;
    end
end
    SER = ser/N;