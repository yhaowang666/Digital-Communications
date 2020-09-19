function [l]=Antenna_mapping(sm1,sm2)
N=length(sm1);
%%s_c=zeros(1,N);s_s=s_c;
for i=1:N
    if sm1(i)==0&&sm2(i)==0
    l(i)=1;
    elseif sm1(i)==0&&sm2(i)==1
        l(i)=2;
            elseif sm1(i)==1&&sm2(i)==0
        l(i)=3;
            elseif sm1(i)==1&&sm2(i)==1
        l(i)=4;
    end
end
    