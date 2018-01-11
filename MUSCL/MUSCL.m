function [ Q_l,Q_r ] = MUSCL( Q )
% MUSCL-TVD
%   Q: original Q
%   Q_MUSCL: used for computing F

omega=1/3; 
omega2=1.5;% 1<=omega2<=(3-omega)/(1-omega)
N=length(Q);
Q_l=zeros(N-1,3);
Q_r=zeros(N-1,3);

% computing Q_l
for i=2:N-1
    Q_l(i,:)=Q(i,:)+0.25*(...
        (1-omega)*LimiterMinmod(Q(i,:)-Q(i-1,:),omega2*(Q(i+1,:)-Q(i,:)))...
        +(1+omega)*LimiterMinmod(Q(i+1,:)-Q(i,:),omega2*(Q(i,:)-Q(i-1,:)))...
        );    
end
Q_l(1,:)=Q(1,:);

% computing Q_r
for i=1:N-2
    Q_r(i,:)=Q(i+1,:)-0.25*(...
        (1-omega)*LimiterMinmod(Q(i+2,:)-Q(i+1,:),omega2*(Q(i+1,:)-Q(i,:)))...
        +(1+omega)*LimiterMinmod(Q(i+1,:)-Q(i,:),omega2*(Q(i+2,:)-Q(i+1,:)))...
        );
end
Q_r(N-1,:)=Q(N,:);
    
end


