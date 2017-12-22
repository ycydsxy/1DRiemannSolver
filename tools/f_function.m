function [ f ] = f_function( p_star,p,rho )
% f(p_star) function of Godunov Scheme
%   p_star: static pressure after the wave
%   p: static pressure before the wave
%   rho: density before the wave
%   f: value of f(p_star)

global gamma R;

T=p/R/rho;% ideal gas
c=sqrt(gamma*R*T);

if p_star>p
    f=(p_star-p)/rho/c/sqrt((gamma+1)/2/gamma*p_star/p+(gamma-1)/2/gamma);
else
    f=2*c/(gamma-1)*((p_star/p).^((gamma-1)/2/gamma)-1);
end

end

