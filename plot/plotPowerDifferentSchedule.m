clc;
clear;
%N0 = 1e-14, tau = 1000~5000
%{ 
BB =         [1.5318e-05 1.3488e-05 1.2929e-05 1.2658e-05 1.2499e-05];
maxSNR =     [1.6811e-05 1.4850e-05 1.4250e-05 1.3959e-05 1.3787e-05];
maxSumRate = [3.3757e-05 3.1248e-05 3.0454e-05 3.0065e-05 2.9834e-05];
greedyPhy =  [1.5268e-05 1.3466e-05 1.2915e-05 1.2648e-05 1.2491e-05];

tier2Time = [10 20 30 40 50]; %s
%}

%N0 = 1e-13, tau = 100~500
BB =         [0.0972 0.0705 0.0328 0.0276 0.0428];
maxSNR =     [0.9084 0.0044 0.0049 ];
maxSumRate = [];
greedyPhy =  [];

tier2Time = [1 2 3 4 5]; %s

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