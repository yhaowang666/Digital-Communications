function b = hamming_encoding(a)
Q = [1 1 1;1 1 0; 1 0 1;0 1 1 ];
G = [eye(4),Q];
N = length(a);
M = zeros(4,N/4);
for i = 1:N
    M(i) = a(i); %matlabĬ���������
end
M = M.';  %M�����ת��
A = M * G;
b = zeros(1, N/4*7);
m = 1;
for i = 1:size(A,1)  %A������
    for j = 1:size(A,2) %A������
        b(m) = A(i,j);
        m = m+1;
    end
end
b = mod(b,2); %ģ2����
end



    