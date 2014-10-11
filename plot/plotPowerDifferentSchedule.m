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

%N0 = 1e-16, tau = 50~100
BB =        fliplr( [5.1921 3.6745 3.5853 3.5729 3.5679] );
maxSNR =    fliplr( [19.6194 9.3930 8.9801 8.9237 8.9013] );
maxEntropy = fliplr( [inf 8.9037 8.1554 8.0616 8.0248] );
greedyPhy = fliplr( [6.8227 3.8674 3.7426 3.7255 3.7187] );

tier2Time = fliplr( [1 10 40 70 100] ); %ms

%plot(xaxis, baselineClusterTier1Power, 'x-.','LineWidth',2,'Color','b','DisplayName',str_baseline ); hold on;
plot(tier2Time,maxSNR,'^-.','LineWidth',2,'DisplayName','Max SNR','Color','g','MarkerSize',10); hold on;
plot(tier2Time,maxEntropy,'*-.','LineWidth',2,'DisplayName','Max entropy','Color','b','MarkerSize',10); hold on;
plot(tier2Time,greedyPhy,'o-.','LineWidth',2,'DisplayName','Greedy physical','Color','c','MarkerSize',10);
plot(tier2Time,BB,'x-.','LineWidth',2,'DisplayName','Branch and bound','Color','r','MarkerSize',10); hold on;
%marker1 = scatter(xMarker,yMarker,'ro');
%set(marker1, 'sizedata', 50);
%set(get(get(marker1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,'XDir','reverse');
ylabel('Tier-2 total power (Watt)');
xlabel('Tier-2 time slot duration (ms)');
legend('show','location','best');
%axis([0 5 -inf inf]);
grid on;