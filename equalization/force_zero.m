function [ c ] = force_zero( h,N )
%%��������ʵ��
%h����һ���Ķྶ�ŵ�ϵ��h
%N����������ͷ��2N+1
%c�����������������ϵ��C
H=length(h);
MID=find(h==1);  %�ҵ�h��ʱ��ԭ��
%��ԭ�������ֵ�������ȣ�����ʹ֮��ȣ�
if(MID-1<H-MID)
    for i=1:(H-MID)-(MID-1)
        h=[0,h];
    end
else
    for i=1:(MID-1)-(H-MID)
        h=[h,0];
    end
end
l=max(MID-1,H-MID);

%���ݸ�����ͷ��ȷ���������x
%x=[h(-2N) h(-2N+1) ... h(0) ... h(2N-1) h(2N)]
x=zeros(1,4*N+1);
if 2*N>=l
    x(2*N+1-l:2*N+1+l)=h;
else
    x=h(MID-2*N:MID+2*N);
end
%����x������󷽳�ϵ��X
%X=[ x(0)     x(-1)    ...    x(-2N) ;
%    x(1)     x(0)     ...    x(-2N+1);
%    ..............................
X=[];
for i=1:2*N+1
    X=[X;fliplr(x(i:2*N+i))]; %fliplr����������ת����
end
%dΪdelta����
d=zeros(2*N+1,1);
d(N+1)=1;
c=X^(-1)*d;
 
end