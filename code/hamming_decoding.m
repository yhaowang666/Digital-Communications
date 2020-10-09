function S = hamming_decoding(a)
Q = [1 1 1;1 1 0; 1 0 1;0 1 1 ];
H = [Q.',eye(3)];
N = length(a);
B = zeros(7,N/7);
for i = 1:N
    B(i) = a(i); %matlabĬ���������
end
B = B.';  %M�����ת��
S = B * (H.');
S = mod(S ,2)
