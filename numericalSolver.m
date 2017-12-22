function result = numericalSolver(S,delta_x,CFL,t_max,schemeType)
% Solver Using Different Schemes
%   S: [rho u p]
%   delta_x: grid spacing
%   CFL: CFL number
%   t_max: total time
%   schemeType: 'Godunov' or 'Roe' etc.
%   result: [rho u p]

global gamma;

N=length(S);
Q=S2Q(S);
F=zeros(N-1,3);

current_time=0;
flag=1;

while(flag)
    delta_t=CFL*delta_x/max(abs(S(:,2))+(sqrt(gamma*S(:,3)./S(:,1))));% timestep
    
    if current_time+delta_t>t_max % stoping criteria
        delta_t=t_max-current_time;
        flag=0;
    end
    current_time=current_time+delta_t;
    
    for i=1:N-1 % computing F_i+1/2
        switch schemeType
            case 'Godunov'
                F(i,:)=godunovScheme(Q(i,:),Q(i+1,:),delta_t);
            case 'Roe'
                
            otherwise
                error('The Scheme You Choose ISNOT Supported!')
        end
    end
    
    for i=2:N-1 % computing Q_i^(n+1) using FISRT ORDER time scheme
        Q(i,:)=Q(i,:)-delta_t/delta_x*(F(i,:)-F(i-1,:));
    end
    
    Q(1,:)=Q(2,:);% fake boundary
    Q(N,:)=Q(N-1,:);% fake boundary
    
    S=Q2S(Q);
end

result=[S(:,1),S(:,2),S(:,3)];
