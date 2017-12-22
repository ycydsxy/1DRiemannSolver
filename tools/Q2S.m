function [ S ] = Q2S( Q )
% transform Q to S
%   Q: [rho rho*u rho*E]
%   S: [rho u p]

global gamma;

s=Q(:,1);
t=Q(:,2);
m=Q(:,3);

rho=s;
u=t./s;
p=(gamma-1).*m;

S=[rho,u,p];

end

