function SER = symbol_error_QPSK(x,y)
N = length(x)/2;
ser = 0;
for i = 1:N
    if x(2*i-1)~= y(2*i-1) || x(2*i)~= y(2*i)
        ser =ser+1;
    end
end
    SER = ser/N;
        