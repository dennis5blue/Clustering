clc;
clear;

recordObj = [Inf, 9.25532013533858e+18, 125477530.908545, 7614.41685514409, 7614.41546707591, 7614.37393262673, 7538.78845805876, ...
    7538.77401184487, 7538.76854462566, 7538.76307764652, 7538.76307566406, 7526.44827941759, 7526.42339595997, 7469.51831000844, ...
    7469.25918401901, 7469.02610584284, 7421.41737887043, 7411.27840970890, 7341.36354012558, 7339.77227099894, 7339.52376419583, ...
    7339.52181990520, 7334.96138721149, 7333.47837047896, 7332.40187353779, 7306.01587431959, 7306.00271168915, 7305.97615702989; ...
    0, 1, 2, 3, 67, 79, 4481, 8712, 49283, 152843, 164363, 2131007, 2359967, 4391999, 4431890, 4858368, 8521768, 8521985, 8528002, ...
    8798272, 9244805, 9245316, 9249407, 10494725, 12849219, 17175693, 17728652, 17859081;];

recordObjNoBound = [Inf, 1.74146052017455e+18, 125477530.908545, 90298593.2585724, 28697202.0510833, 7614.40162031440, ...
    7614.40024043796, 7614.37947198728, 7538.77537496337, 7538.77537496337, 7538.77537496337; ...
    0, 1, 2, 3, 4, 5, 7, 16, 39, 641 5001;];
%recordObj = load(['./data/BBconverge.mat']);
%recordObjNoBound = load(['./data/BBNoBoundconverge.mat']);
startPoint = 3;

numOfIteration = 1000;
iteration = 1:numOfIteration;
vecBBPayoff = recordObj(1,startPoint:length(recordObj(1,:)));
vecBBIter = recordObj(2,startPoint:length(recordObj(2,:)));
vecBBPlot = zeros(1,numOfIteration);

vecBBPayoff2 = recordObjNoBound(1,startPoint:length(recordObjNoBound(1,:)));
vecBBIter2 = recordObjNoBound(2,startPoint:length(recordObjNoBound(2,:)));
vecBBPlot2 = zeros(1,numOfIteration);

minPayoff = vecBBPayoff(length(vecBBPayoff));

xMarker = vecBBIter;
yMarker = vecBBPayoff/minPayoff;
xMarker2 = vecBBIter2;
yMarker2 = vecBBPayoff2/minPayoff;

for i=1:numOfIteration
    temp = find(vecBBIter>=i);
    payoff = vecBBPayoff( temp(1) );
    vecBBPlot(i) = payoff/minPayoff;
    
    temp2 = find(vecBBIter2>=i);
    payoff2 = vecBBPayoff2( temp2(1) );
    vecBBPlot2(i) = payoff2/minPayoff;
end
%DC_colorMap=jet(16);

semilogx(iteration,vecBBPlot(1:numOfIteration),'LineWidth',3,'DisplayName','Branch and bound','Color','r','MarkerSize',10); hold on;
marker1 = scatter(xMarker,yMarker,'ro');
set(marker1, 'sizedata', 50);
set(get(get(marker1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

semilogx(iteration,vecBBPlot2(1:numOfIteration),'LineWidth',3,'DisplayName','Full search','Color','b','MarkerSize',10);
marker2 = scatter(xMarker2,yMarker2,'bo');
set(marker2, 'sizedata', 50);
set(get(get(marker2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

ylabel('Objective Value (%)');
xlabel('Iteration');
legend('show');
axis([0 numOfIteration -inf inf]);
grid on;