clc;
clear;

totalNodes=20;
N = 4;
totalTime = 100.0; % 100ms

BB =         [2.0140 0.5161 0.0103 1.0275];
maxSNR =     [7.4914 1.2771 0.0601 0.0727];
maxEntropy = [2.6090 0.9546 4.3686 0.926];
greedyPhy =  [2.2870 1.3287 0.0926 0.0103];

%BB = sort(BB,'descend');
%maxSNR = sort(maxSNR,'descend');
%maxEntropy = sort(maxEntropy,'descend');
%greedyPhy = sort(greedyPhy,'descend');

slot = [1 2 3]; 

tier2Power = [];
for i=1:length(slot)
    tier2Power = [tier2Power; maxSNR(i) maxEntropy(i) greedyPhy(i) BB(i)];
end

subplot(2,1,1);
Obar = bar(slot,tier2Power,'hist');
set(Obar(1),'FaceColor',[1 0.35 0.35]);
set(Obar(2),'FaceColor',[0 1 0.7]);
set(Obar(3),'FaceColor',[1 1 0.7]);
set(Obar(4),'FaceColor',[0 1 1]);

set(gca,'xticklabel');
xlabel(gca,'Slot index');
ylabel(gca,'Tier-2 Power (Watt)');
%title('|S|=24, |H|=4, T+N*\tau = 100(s)');
legend('Max SNR','Max entropy','Greedy physical','Branch and bound','location','Best');
grid on;
hold off;
