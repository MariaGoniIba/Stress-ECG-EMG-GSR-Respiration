function [Nr,Ni,E]=recurrence_analysis(x)

X=[x(1:end-4) x(2:end-3) x(3:end-2) x(4:end-1) x(5:end)];
D=dist(X,X');
I=D<0.05*mean(mean(D));

%Delete diagonal
I=I.*(1-eye(size(I,1)));

N=zeros(1,4);
N(1)=sum(sum(I)); %Number of total points

n=1;
while(N(n)~=0)
    n=n+1;
    I=I(2:end,1:end-1).*I(1:end-1,2:end);
    N(n)=sum(sum(I));
end

Nr=N(1)/size(D,1)^2;
Ni=N(3)/size(D,1)^2;
P=N(1:end-1)-N(2:end);
P=P/sum(P);
E=-sum(P.*log2(P+eps));

