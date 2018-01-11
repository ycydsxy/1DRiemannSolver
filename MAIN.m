clc;
clear;
path(path,'.\tools');
path(path,'.\schemes');
path(path,'.\MUSCL');

global gamma R;
gamma=1.4;
R=287;

% NOTICE:
% S: BASIC VARIABLES,[rho u p]
% Q: CONSERVATIVE VARIABLES,[rho rho*u rho*E]
% F: FLUX VARIABLES,[rho*u u^2+p rho*u*H]

N=100;% grid number
CFL=0.5;% CFL number
x=linspace(0,1,N)';% grid
x1=linspace(0,1,1000)';% grid for exact solution
delta_x=(x(end)-x(1))/(N-1);% grid spacing
x_0=0.5;% discontinuity surface location
n_0=ceil((x_0-x(1))/delta_x);% index of discontinuity surface
totalTimes=[0.2,0.15,0.012,0.035,0.035];

exactSolutions=cell(1,5);
numericalSolutions=cell(5,5);

test_num=[2];

for i=test_num
    S=initialize(N,n_0,i);% initial field of TESTi
    t_max=totalTimes(i);

    exactSolutions{i}=exactRiemannSovler(S(1,:),S(end,:),x1,t_max,x_0);
    %numericalSolutions{1,i}=numericalSolver(S,delta_x,CFL,t_max,'Godunov',0);% Godunov Scheme
    numericalSolutions{2,i}=numericalSolver(S,delta_x,CFL,t_max,'Roe',0);% Roe Scheme
    numericalSolutions{3,i}=numericalSolver(S,delta_x,CFL,t_max,'AUSM',0);% AUSM Scheme
    %numericalSolutions{4,i}=numericalSolver(S,delta_x,CFL,t_max,'Roe',1);% Roe Scheme + MUSCL
    %numericalSolutions{5,i}=numericalSolver(S,delta_x,CFL,t_max,'AUSM',1);% AUSM Scheme + MUSCL
end






