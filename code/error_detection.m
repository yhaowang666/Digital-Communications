function y = error_detection(a,S)
for i = 1 : size(S,1)
    if S(i,:) == [1 1 1]
        a((i-1)*7+1) = mod(a((i-1)*7+1)+1 , 2);
    elseif S(i,:) == [1 1 0]
        a((i-1)*7+2) = mod(a((i-1)*7+2)+1 , 2);
    elseif S(i,:) == [1 0 1]
        a((i-1)*7+3) = mod(a((i-1)*7+3)+1 , 2);
    elseif S(i,:) == [0 1 1]
        a((i-1)*7+4) = mod(a((i-1)*7+4)+1 , 2);
    elseif S(i,:) == [1 0 0]
        a((i-1)*7+5) = mod(a((i-1)*7+5)+1 , 2);
    elseif S(i,:) == [0 1 0]
        a((i-1)*7+6) = mod(a((i-1)*7+6)+1 , 2);    
    elseif S(i,:) == [0 0 1]
        a((i-1)*7+7) = mod(a((i-1)*7+7)+1 , 2);
    end
end
y = a;

        