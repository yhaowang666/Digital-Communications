function y = hamming_decoding(a)
Q = [1 1 1;1 1 0; 1 0 1;0 1 1 ];
H = [Q.',eye(3)];
N = length(a);
B = zeros(7,N/7);
for i = 1:N
    B(i) = a(i); %matlab默认纵向填充
end
B = B.';  %M矩阵的转置
S = B * (H.');
S = mod(S ,2); %模2运算
y1 = error_detection(a,S);  %纠错
for i = 0:length(y1)/7-1
    y(i*4+1:i*4+4) = y1(i*7+1:i*7+4);  %取码元（7bit）中前4bit为译码输出
end

        
        
 

