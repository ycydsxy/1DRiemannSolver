function result = numericalSolver(S,delta_x,CFL,t_max,schemeType,M_flag)
% Solver Using Different Schemes
%   S: [rho u p]
%   delta_x: grid spacing
%   CFL: CFL number
%   t_max: total time
%   schemeType: 'Godunov' or 'Roe' etc.
%   M_flag:use MUSCL or not 
%   result: [rho u p]

tic;

global gamma;

N=length(S);
Q=S2Q(S);
F=zeros(N-1,3);

current_time=0;
steps=0;
flag=1;
maxSteps=1e5;

while(flag)
    steps=steps+1;
    delta_t=CFL*delta_x/max(abs(S(:,2))+(sqrt(gamma*S(:,3)./S(:,1))));% timestep
    
    if steps>maxSteps % stoping criteria
        disp('WARNING: maxSteps reached!');
        break;
    end
    
    if current_time+delta_t>t_max % stoping criteria
        delta_t=t_max-current_time;
        flag=0;
    end
    current_time=current_time+delta_t;
    
    if M_flag
       [Q_l,Q_r]=MUSCL( Q );
    else
        Q_l=Q(1:N-1,:);
        Q_r=Q(2:N,:);
    end
    
    
    for i=1:N-1 % computing F_i+1/2
        switch schemeType
            case 'Godunov'
                F(i,:)=godunovScheme(Q_l(i,:),Q_r(i,:),delta_t);
            case 'Roe'
                F(i,:)=roeScheme(Q_l(i,:),Q_r(i,:));
            otherwise
                error('The Scheme You Choose ISNOT Supported!')
        end
    end
    
    for i=2:N-1 % computing Q_i^(n+1) using FISRT ORDER time scheme
        Q(i,:)=Q(i,:)-delta_t/delta_x*(F(i,:)-F(i-1,:));
    end
    
   % no boundary conditions needed
    
    S=Q2S(Q);    
end
toc
result=[S(:,1),S(:,2),S(:,3)];
