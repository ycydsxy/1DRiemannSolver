%% clear null cells
numericalSolutions(cellfun(@isempty,numericalSolutions))=[];
%% ploting
options={'o','s','d','^','v'};
ylabels={'\rho','u','p'};

for i=1:5
    figure
    result=exactSolutions{1,i};
    for j=1:3
        subplot(3,1,j);
        plot(x1,result(:,j));
        if j==1
            title(sprintf('TEST%d t=%.3f,CFL=%.2f,N=%d',i,t_max,CFL,N));
        end
        
        for k=1:size(numericalSolutions)
            result_numerical=numericalSolutions{k,i};
            hold on
            plot(x,result_numerical(:,j),options{k});
        end
        legend('Exact Solution','Godunov Scheme','Roe','Roe\_MUSCLE','AUSM');
        xlabel('x');
        ylabel(ylabels{j});
    end
end