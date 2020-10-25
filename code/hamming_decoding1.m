function y = hamming_decoding1(x)
A = [0 0 0 0 0 0 0;
     0 0 0 1 0 1 1;
     0 0 1 0 1 0 1;
     0 0 1 1 1 1 0;
     0 1 0 0 1 1 0;
     0 1 0 1 1 0 1;
     0 1 1 0 0 1 1;
     0 1 1 1 0 0 0;
     1 0 0 0 1 1 1;
     1 0 0 1 1 0 0;
     1 0 1 0 0 1 0;
     1 0 1 1 0 0 1;
     1 1 0 0 0 0 1;
     1 1 0 1 0 1 0;
     1 1 1 0 1 0 0;
     1 1 1 1 1 1 1];
     B = -sqrt(2).*A + sqrt(2)/2;
     for i = 0 : length(x)/7-1
         m = x(7*i+1:7*i+7);
         d = sum((m - B(1,:)).^2);
         index = 1;
         for j = 2 : 16
             n = sum((m - B(j,:)).^2);
             if (d > n)
                d = n;
                index = j;
             end
         end
         y(4*i+1:4*i+4) = A(index,1:4);
     end
     
            
             
         