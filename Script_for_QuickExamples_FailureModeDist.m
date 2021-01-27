%Wiebull Plots

MTTF = [6000 3500 2500]
shape = [5.6 5.0 2.9]
stoptime = 2700;
t = 1:3e3;
x=zeros(length(t), length(MTTF))
for d = 1:length(MTTF)
    x(:,d) = wblpdf(t,MTTF(d),shape(d));
    nm{d} = sprintf('Failure Mode %i',d);
end
figure(1)
subplot(212)
plot(t,x,'-')
hold on
xline(stoptime);
title('Individual Failure Mode Distributions')
xlabel('Minutes')
ylabel('Probability of Failure')
legend(nm{:},'Simulation Stop','location','northwest')
hold off
subplot(211)
plot(t,mean(x,2))
xline(stoptime);
hold on;
title('Machine Failure Distribution')
xlabel('Minutes')
ylabel('Probability of Failure')
legend('Total Failure Rate','Simulation Stop','location','northwest')
hold off;



figure(2)
plot(t,x,'-',t,mean(x,2),':')
hold on
xline(stoptime);
title('Individual Failure Mode & Total Machine Failure Distributions')
xlabel('Minutes')
ylabel('Failure Distribution')
legend(nm{:},'Simulation Stop','Total Failure Rate','location','northwest')

hold off