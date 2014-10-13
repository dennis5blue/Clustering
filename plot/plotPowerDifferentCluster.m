clc;
clear;

totalNodes=20;
N = 4;
%{
totalTime = 1000.0; % 1s
CSAclusterTier1Power = [0.0057 0.0047 0.0040 0.0035 0.0031];
CSAclusterTier2Power = [3.5656 3.5679 3.5718 3.5795 3.6028];

KmeansclusterTier1Power = [0.0163 0.0136 0.0116 0.0102 0.0090];
KmeansclusterTier2Power = [0.7002 1.0932 3.5536 131.7352 1.2699e+09];

DMCPclusterTier1Power = [0.5892 0.0295 0.0252 0.0221 0.0196];
DMCPclusterTier2Power = [4.5227 9.6166 83.9697 1.6979e+05 5.4372e+18];

xaxis = [500.0 600.0 700.0 800.0 900.0]; % ms
%}

totalTime = 1000.0; % 1s
CSAclusterTier1Power = [0.7107 0.5668 0.4714 0.4304];
CSAclusterTier2Power = [0.4912 0.5900 0.7019 0.7580];

KmeansclusterTier1Power = [2.0483 1.6335 1.3584 1.1627];
KmeansclusterTier2Power = [0.5512 0.8019 1.5314 5.7637];

DMCPclusterTier1Power = [4.4447 3.5447 2.9478 2.5230];
DMCPclusterTier2Power = [2.9987 4.5227 9.6166 83.9697];

xaxis = [400.0 500.0 600.0 700.0]; % ms

subplot(2,1,1);
str_DMCP = sprintf('DMCP');
semilogy(xaxis, DMCPclusterTier1Power, '*-','LineWidth',2,'Color','g','DisplayName',str_DMCP ); hold on;
str_baseline = sprintf('Kmeans');
semilogy(xaxis, KmeansclusterTier1Power, 'x-.','LineWidth',2,'Color','b','DisplayName',str_baseline ); hold on;
str_CSA = sprintf('Proposed');
semilogy(xaxis, CSAclusterTier1Power, 'o-','LineWidth',2,'Color','r','DisplayName',str_CSA ); hold on;

xlabel('Tier-1 resource (ms)');
ylabel('Tier-1 Power (Watt)');
axis([-inf inf 0.4 5.0]);
legend('show','location','NorthEast');
grid on;

subplot(2,1,2);
str_DMCP = sprintf('DMCP');
semilogy(xaxis, DMCPclusterTier2Power, '*-','LineWidth',2,'Color', 'g', 'DisplayName',str_DMCP); hold on;
str_baseline = sprintf('Kmeans');
semilogy(xaxis, KmeansclusterTier2Power, 'x-.','LineWidth',2,'Color', 'b', 'DisplayName',str_baseline); hold on;
str_CSA = sprintf('Proposed');
semilogy(xaxis, CSAclusterTier2Power, 'o-','LineWidth',2,'Color', 'r', 'DisplayName',str_CSA); hold on;

xlabel('Tier-1 resource (ms)');
ylabel('Tier-2 Power (Watt)');
axis([-inf inf 0.4 inf]);
%title('|S|=24, |H|=4, T+N*\tau = 100(s)');
legend('show','location','NorthWest');
grid on;
hold off;
