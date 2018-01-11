%% clear null cells
numericalResults=numericalSolutions(2:3,:);
%% ploting
options={'-.or','-.sg','-.db','-.^c','-.vm'};
ylabels={'\rho','u','p'};

for i=test_num
    result=exactSolutions{1,i};
    figure;
    for j=1:3
        subplot(3,1,j);
        plot(x1,result(:,j),'k');
        if j==1
            title(sprintf('TEST%d t=%.3f,CFL=%.2f,N=%d',i,t_max,CFL,N));
        end
        for k=1:size(numericalResults)
            result_numerical=numericalResults{k,i};
            hold on
            plot(x,result_numerical(:,j),options{k},'MarkerSize',4);
        end
        legend('Exact Solution','Roe','AUSM','Roe\_MUSCLE','AUSM\_MUSCLE');
        %legend('Exact Solution','Godunov Scheme','Roe','AUSM','Roe\_MUSCLE','AUSM\_MUSCLE');
        xlabel('x');
        ylabel(ylabels{j});
    end
end