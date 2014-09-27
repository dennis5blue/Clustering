clc;
clear;
%{
BB = [1.4754e-10 1.3038e-10 1.2512e-10 1.2257e-10 1.2107e-10];
maxSNR = [1.6193e-10 1.4357e-10 1.3792e-10 1.3518e-10 1.3356e-10];
maxSumRate = [3.2578e-10 3.0224e-10 2.9477e-10 2.9110e-10 2.8893e-10];
greedyPhy = []
%}
BB =         [1.5318e-05 1.3488e-05 1.2929e-05 1.2658e-05 1.2499e-05];
maxSNR =     [1.6811e-05 1.4850e-05 1.4250e-05 1.3959e-05 1.3787e-05];
maxSumRate = [3.3757e-05 3.1248e-05 3.0454e-05 3.0065e-05 2.9834e-05];
greedyPhy =  [1.5268e-05 1.3466e-05 1.2915e-05 1.2648e-05 1.2491e-05];

tier2Time = [10 20 30 40 50]; %ms

%plot(xaxis, baselineClusterTier1Power, 'x-.','LineWidth',2,'Color','b','DisplayName',str_baseline ); hold on;

plot(tier2Time,BB,'x-.','LineWidth',2,'DisplayName','Branch and bound','Color','r','MarkerSize',10); hold on;
plot(tier2Time,maxSNR,'^-.','LineWidth',2,'DisplayName','Max SNR','Color','g','MarkerSize',10); hold on;
plot(tier2Time,maxSumRate,'*-.','LineWidth',2,'DisplayName','Max sum rate','Color','b','MarkerSize',10); hold on;
plot(tier2Time,greedyPhy,'o-.','LineWidth',2,'DisplayName','Greedy physical','Color','c','MarkerSize',10);
%marker1 = scatter(xMarker,yMarker,'ro');
%set(marker1, 'sizedata', 50);
%set(get(get(marker1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

ylabel('Tier-2 power (Watt)');
xlabel('Time slot duration (s)');
legend('show','location','best');
%axis([0 5 -inf inf]);
grid on;