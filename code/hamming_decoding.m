function y = hamming_decoding(a)
Q = [1 1 1;1 1 0; 1 0 1;0 1 1 ];
H = [Q.',eye(3)];
N = length(a);
B = zeros(7,N/7);
for i = 1:N
    B(i) = a(i); %matlabĬ���������
end
B = B.';  %M�����ת��
S = B * (H.');
S = mod(S ,2); %ģ2����
y1 = error_detection(a,S);  %����
for i = 0:length(y1)/7-1
    y(i*4+1:i*4+4) = y1(i*7+1:i*7+4);  %ȡ��Ԫ��7bit����ǰ4bitΪ�������
end

        
        
 

