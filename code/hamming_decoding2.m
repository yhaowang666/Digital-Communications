function y = hamming_decoding2(a)
N = length(a);
B = zeros(7,N/7);
for i = 1:N
    B(i) = a(i); %matlab默认纵向填充
end
B = B.';  %M矩阵的转置
y = mhd(B);    %最小汉明距离译码


        
        
 

