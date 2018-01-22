%% clear null cells
numericalResults=numericalSolutions(1:4,:);
%% ploting
options={'-.or','-.sg','-.db','-.^m','-.vc'};
ylabels={'\rho','u','p'};

for i=test_num
    result=exactSolutions{1,i};
    figure;
    for j=1:3
        subplot(3,1,j);
        plot(x1,result(:,j),'k');
        if j==1
            title(sprintf('TEST%d x_0=%0.2f,t=%.3f,CFL=%.2f,N=%d',i,x_0,t_max,CFL,N));
        end
        for k=1:size(numericalResults)
            result_numerical=numericalResults{k,i};
            hold on
            plot(x,result_numerical(:,j),options{k},'MarkerSize',4);
        end
        legend('Exact Solution','Godunov','Roe','AUSM','Roe\_MUSCLE');
        xlabel('x');
        ylabel(ylabels{j});
    end
end