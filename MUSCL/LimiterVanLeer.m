function [ phi ] = LimiterVanLeer( a,b )
% Van Leer Limiter

N=length(a);
N1=length(b);

if N~=N1
    error('Dimension Error!');
end

phi=zeros(1,N);

for i=1:N
    
flag=a(i)*b(i);

if flag<=0
    phi(i)=0;
else
    phi(i)=(flag+abs(flag))/(a(i)+b(i));
end
end
end

