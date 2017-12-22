clc;
clear;
path(path,'.\tools');

global gamma R;
gamma=1.4;
R=287;

% NOTICE:
% S: BASIC VARIABLES,[rho u p]
% Q: CONSERVATIVE VARIABLES,[rho rho*u rho*E]
% F: FLUX VARIABLES,[rho*u u^2+p rho*u*H]

N=20;% grid number
CFL=0.9;% CFL number
x=linspace(-1,1,N)';% grid
x1=linspace(-1,1,2000)';% grid for exact solution
delta_x=2/(N-1);% grid spacing
x_0=0;% discontinuity surface location
n_0=floor((x_0+1)/delta_x+1);% index of discontinuity surface
totalTimes=[0.25,0.15,0.012,0.035,0.012];

exactSolutions=cell(1,5);
numericalSolutions=cell(5,5);

for i=1:5
    S=initialize(N,n_0,i);% initial field of TESTi
    t_max=totalTimes(i);

    exactSolutions{i}=exactRiemannSovler(S(1,:),S(end,:),x1,t_max);
    numericalSolutions{1,i}=numericalSolver(S,delta_x,CFL,t_max,'Godunov');
    
end

%ploting
ploting;






