function S = hamming_decoding(a)
Q = [1 1 1;1 1 0; 1 0 1;0 1 1 ];
H = [Q.',eye(3)];
N = length(a);
B = zeros(7,N/7);
for i = 1:N
    B(i) = a(i); %matlab默认纵向填充
end
B = B.';  %M矩阵的转置
S = B * (H.');
S = mod(S ,2)
