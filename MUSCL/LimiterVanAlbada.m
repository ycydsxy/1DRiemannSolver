function [ phi ] = LimiterVanAlbada( a,b )
% VanAlbada Limiter

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
    phi(i)=flag*(a(i)+b(i))/(a(i)^2+b(i)^2);
end
end
end

