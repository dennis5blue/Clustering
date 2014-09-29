clc;
clear;

totalNodes=20;
N = 4;
totalTime = 10000.0; % 10s

CSAclusterTier1Power = [0.0208 0.0160 0.0130 0.0110 0.0095];
CSAclusterTier2Power = [1.4557e-05 1.5318e-05 1.6678e-05 1.9793e-05 3.1634e-05];
CSAtotalPower = CSAclusterTier1Power + CSAclusterTier2Power;

KmeansclusterTier1Power = [0.0589 0.0456 0.0371 0.0313 0.0271];
KmeansclusterTier2Power = [0.2069 0.2129 0.2235 0.2467 0.3378];
KmeansTotalPower = KmeansclusterTier1Power + KmeansclusterTier2Power;

DMCPclusterTier1Power = [0.1215 0.0956 0.0787 0.0669 0.0581];
DMCPclusterTier2Power = [0.8118 0.8425 0.8968 1.0191 1.5411];
DMCPTotalPower = DMCPclusterTier1Power + DMCPclusterTier2Power;

xaxis = [5.0 6.0 7.0 8.0 9.0 ]; % sec

str_DMCP = sprintf('DMCP');
plot(xaxis, DMCPTotalPower, '*-','LineWidth',2,'Color','g','DisplayName',str_DMCP ); hold on;
str_baseline = sprintf('Kmeans clustering algorithm');
plot(xaxis, KmeansTotalPower, 'x-.','LineWidth',2,'Color','b','DisplayName',str_baseline ); hold on;
str_CSA = sprintf('Proposed clustering algorithm');
plot(xaxis, CSAtotalPower, 'o-','LineWidth',2,'Color','r','DisplayName',str_CSA ); hold on;

xlabel('Tier 1 resource (s)');
ylabel('Total Power (Watt)');
title('|S|=20, |H|=4, T+N*\tau = 10(s)');
legend('show','location','NorthWest');
grid on;

hold off;
