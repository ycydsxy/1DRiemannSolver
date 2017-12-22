function [ S ] = initialize( N,n_0,index )
% Initialize fuild field for different TESTs
%   N: grid number
%   n_0: index of discontinuity surface
%   index: TEST index
%   S: [rho u p]


switch index
    case 1
        S_l=[1,0,1];
        S_r=[0.125,0,0.1];
    case 2
        S_l=[1,-2,0.4];
        S_r=[1,2,0.4];
    case 3
        S_l=[1,0,1000];
        S_r=[1,0,0.01];
    case 4
        S_l=[1,0,0.01];
        S_r=[1,0,100];
    case 5
        S_l=[5.99924,19.5975,460.894];
        S_r=[5.99242,-6.19633,46.0950];
end

for i=1:N
    if i<=n_0
        S(i,:)=S_l;
    else
        S(i,:)=S_r;
    end
end
end

