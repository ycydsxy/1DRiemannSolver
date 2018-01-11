function F = godunovScheme(Q_l,Q_r,t)
% Computing F_i+1/2 Using Godunov Scheme
%   Q_l: [rho rho*u rho*E] on the left side
%   Q_r: [rho rho*u rho*E] on the right side
%   t: time
%   F: [rho*u u^2+p rho*u*H] 

global gamma R;

S_l=Q2S(Q_l);
S_r=Q2S(Q_r);

rho_1=S_l(:,1);
u_1=S_l(:,2);
p_1=S_l(:,3);
rho_2=S_r(:,1);
u_2=S_r(:,2);
p_2=S_r(:,3);

F_0=f_function(0,p_1,rho_1)+f_function(0,p_2,rho_2);
F_p1=f_function(p_1,p_1,rho_1)+f_function(p_1,p_2,rho_2);
F_p2=f_function(p_2,p_1,rho_1)+f_function(p_2,p_2,rho_2);
delta_u=u_1-u_2;

% determine which condition the problem belongs to
if delta_u<=F_0
    condition=5;
elseif delta_u>F_0 && delta_u<=min(F_p1,F_p2)
    condition=4;
elseif delta_u>max(F_p1,F_p2)
    condition=1;
else
    if p_1>p_2
        condition=2;
    else
        condition=3;
    end
end

%solving equation using Newton Method
fun=@(p_star)f_function(p_star,p_1,rho_1)+f_function(p_star,p_2,rho_2)-delta_u;
[p_star,fval]=fsolve(fun,0.5*(p_1+p_2),optimset('Display','off','TolFun',1e-10));
u_star=0.5*(u_1+u_2+f_function(p_star,p_2,rho_2)-f_function(p_star,p_1,rho_1));
% fprintf('fval = %.4e\n',fval);

% constant parameters
T_1=p_1/R/rho_1;% ideal gas
c_1=sqrt(gamma*R*T_1);
A_1=rho_1*c_1*sqrt((gamma+1)/2/gamma*p_star/p_1+(gamma-1)/2/gamma);
T_2=p_2/R/rho_2;% ideal gas
c_2=sqrt(gamma*R*T_2);
A_2=rho_2*c_2*sqrt((gamma+1)/2/gamma*p_star/p_2+(gamma-1)/2/gamma);

switch condition % computing left wave
    case {1,3}
        rho_1_star= rho_1*A_1/(A_1-rho_1*(u_1-u_star));
        z_1h=u_1-A_1/rho_1;
        z_1t=z_1h;
    case {2,4}
        c_1_star=c_1+(gamma-1)*(u_1-u_star)/2;
        rho_1_star=gamma*p_star/c_1_star^2;
        z_1h=u_1-c_1;
        z_1t=u_star-c_1_star;
    case{5}
        c_1_star=c_1+(gamma-1)*(u_1-u_star)/2;
        rho_1_star=gamma*p_star/c_1_star^2;
        z_1h=u_1-c_1;
        z_1t=u_1-2/(gamma-1)/c_1;
end

switch condition % computing right wave
    case {1,2}
        rho_2_star= rho_2*A_2/(A_2-rho_2*(u_star-u_2));
        z_2h=u_2+A_2/rho_2;
        z_2t=z_2h;
    case {3,4}
        c_2_star=c_2+(gamma-1)*(u_star-u_2)/2;
        rho_2_star=gamma*p_star/c_2_star^2;
        z_2h=u_2+c_2;
        z_2t=u_star+c_2_star;
    case{5}
        c_2_star=c_2+(gamma-1)*(u_star-u_2)/2;
        rho_2_star=gamma*p_star/c_2_star^2;
        z_2h=u_2+c_2;
        z_2t=u_2+2/(gamma-1)/c_2;
end

%fprintf('[p* u* rho1* rho2*]=[%.5f %.5f %.5f %.5f]\n',p_star,u_star,rho_1_star,rho_2_star);

x=0;

if x<z_1h*t
    u=u_1;
    p=p_1;
    rho=rho_1;
elseif x>=z_1h*t && x<z_1t*t
    c_i=(gamma-1)/(gamma+1)*(u_1-x/t)+2/(gamma+1)*c_1;
    u=x/t+c_i;
    p=p_1*(c_i/c_1)^(2*gamma/(gamma-1));
    rho=gamma*p/c_i^2;
elseif x>=z_1t*t && x<z_2t*t
    u=u_star;
    p=p_star;
    if x<u_star*t
        rho=rho_1_star;
    else
        rho=rho_2_star;
    end
elseif x>=z_2t*t && x<z_2h*t
    c_i=(gamma-1)/(gamma+1)*(x/t-u_2)+2/(gamma+1)*c_2;
    u=x/t-c_i;
    p=p_2*(c_i/c_2)^(2*gamma/(gamma-1));
    rho=gamma*p/c_i^2;
elseif x>=z_2h*t
    u=u_2;
    p=p_2;
    rho=rho_2;
end

F=S2F([ rho,u,p ]);
end

