function y = hamming_decoding2(a)
N = length(a);
B = zeros(7,N/7);
for i = 1:N
    B(i) = a(i); %matlabĬ���������
end
B = B.';  %M�����ת��
y = mhd(B);    %��С������������


        
        
 

