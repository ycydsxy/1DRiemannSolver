function [ S ] = F2S( F )
% transform F to S
%   F: [rho*u rho*u^2+p rho*u*H]
%   S: [rho u p]

global gamma;

s=F(:,1);
t=F(:,2);
m=F(:,3);

p=(t+sqrt(t.^2-4*(gamma-1)/gamma.*m.*s))/2;
u=m*(gamma-1)/gamma./p;
rho=s./u;

S=[rho,u,p];

end

