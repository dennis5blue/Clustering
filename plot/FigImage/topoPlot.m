%   Write in script form first
%   Detailed explanation goes here
%   Mapfile and structure also include path
%   e.g.   showuUL2DClusterStructure( 'mapFile\mapFile_uni_exp_40\mapFile_uni_exp_40_1.txt', 'clusterStructure/Nodes40/DLMatlab_40H15C1.txt',1 )
%   ctrl = 0 -> no line connect with member 
%   ctrl = 1 -> line connect with member
%   Note: almost the same as DL. Only parameter changes when input Power
clear all;close all;
figure;
%dpath = ['data/4b2DC_m13.txt';'data/2i2MC_m9.txt ';'data/SKMDC_m11.txt';'data/SKMMC_m11.txt'];
%strtitle=['Two-Tier DC';'Two-Tier MC';'K-Means  DC';'K-Means  MC'];


CSPath = ['./data/CS_test_CSA_20_v4.out'];
SchedulePath = ['./data/imageSolution_20cam_v4.out'];
mapString = ['./data/20cam_r500_map.out'];
matSolution = dlmread(SchedulePath);
vecDirection = dlmread('./data/20cam_r500_dir.out');

[totalNodes radius x y Gij Gib]= parse_Map_v1(mapString);
radius = 500;
[ maxChNum powerMax C2WT payoffs SAFac TimeMs1st TimeMs2nd EnergyJ1st ...
    EnergyJ2nd clusterStru headList TxPower  ] = parse_Stru_v2 ...
    (CSPath,totalNodes);

nodeIdx = 1:totalNodes;

for ii=1:4
% 
    hFig = figure(ii);
    set(hFig, 'Position', [1 1 550 500 ]);

    powerBound = 1e-6;

    headName=headList;

    densityScale = 1000;%km square
    density = totalNodes / (pi*(radius/densityScale)^2) ;

    hold on;
    %circle(-10,55,radius/1.7);
    circle(0,0,radius/1.2);

    B = find(sum(clusterStru));
    supSet = setdiff(B, headList);
    Binary_Unsupset=ones(1,totalNodes)-sum(clusterStru);
    index_unsupset=find(Binary_Unsupset);
    
    inactiveSet = intersect(supSet, nodeIdx(matSolution(ii,:) == 0));
    activeSet = nodeIdx(matSolution(ii,:) == 1);
    %circle(-10,55,radius/1.7);
    circle(0,0,radius/1.2);
    clusterSize=zeros(1,maxChNum);

    firstTierX(1)=0;
    firstTiery(1)=0;
    %connect the head and serves nodes
    for i=1:maxChNum
        if(headName(i)==0)continue;end
        firstTierX(2)=x(headName(i));
        firstTiery(2)=y(headName(i));
        hLine=plot(firstTierX,firstTiery,'k-','LineWidth',2.0);
        set(get(get(hLine,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
    
        for j=1:totalNodes
            if clusterStru(i,j)==1
                X(1) = x(headName(i));
                Y(1) = y(headName(i));
                X(2) = x(j);
                Y(2) = y(j);
                theta = vecDirection(j)+pi/2;
                r = 50;
                u = x(j) + r * cos(theta); % convert polar (theta,r) to cartesian
                v = y(j) + r * sin(theta);
                if matSolution(ii,j) == 1
                    ARR = arrow([x(j) y(j)],[u v],'BaseAngle',30,'length',15,'width',1.5,...
                        'EdgeColor',[1 0 0],'FaceColor',[1 0 0]);
                    set( get(get(ARR,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                    hLine = plot(X,Y, ...
                        'Color',[1 0 0], ...
                        'LineWidth',2);
                    set( get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                    clusterSize(i) = clusterSize(i)+1;
                else
                    ARR = arrow([x(j) y(j)],[u v],'BaseAngle',30,'length',15,'width',0.4,...
                        'EdgeColor',[0 0 1],'FaceColor',[0 0 1]);
                    set( get(get(ARR,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                    hLine = plot(X,Y,'k:', ...
                        'Color',[0.001 0.001 0.001], ...
                        'LineWidth',2);
                    set( get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                    clusterSize(i) = clusterSize(i)+1;
                end
            end
        end
    end

    dBToRaw = @(dB) 10.^(dB/10);
    [X1, X2] = meshgrid(x(1:totalNodes));
    [Y1, Y2] = meshgrid(y(1:totalNodes));
    X_sq = (X1-X2).^2;
    Y_sq = (Y1-Y2).^2;
    D2 = X_sq + Y_sq;
    Gij = dBToRaw(-(131.1 + 42.81*log10( (D2.^0.5)/1000 )));
    
    %Compute the SINR below
    [Rate SINR Inter]  = ...
        calSINR_4All( maxChNum, totalNodes, clusterStru,  clusterSize, ...
        headName, TxPower,  Gij  );
    
    %Add text
    for i=1:maxChNum
        if(headName(i)==0) continue;end
        for j=1:totalNodes
            if (headName(i)==j && clusterStru(i,j)==1)
                %str = sprintf('H %d', headName(i));
                %text(x(j)+10,y(j)-10,str,'FontSize',12);
                %str = sprintf('|C%d|=%d',i, clusterSize(i)-1);
                %text(x(j)+10,y(j)-10,str,'FontSize',10);
            elseif clusterStru(i,j)==1
                str = sprintf('%.5f',SINR(j));
                %text(x(j),y(j),str,'FontSize',12);
            end
        end
    end

plot(x(index_unsupset),y(index_unsupset),...
    'LineWidth',3, ...
    'MarkerFaceColor',[0.2 1 0], ...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',10,...
    'Marker','X',...
    'LineStyle','none', ...
    'DisplayName','Unselected');

plot(x(inactiveSet),y(inactiveSet), ...
    'MarkerFaceColor', [0 0 1], ...
    'MarkerEdgeColor',[0 0 0], ...
    'Marker', 'o',...
    'LineStyle','none', ...
    'MarkerSize', 10, ...
    'DisplayName','Inactive');


plot(x(activeSet),y(activeSet), ...
    'MarkerFaceColor', [1 0 0], ...
    'MarkerEdgeColor',[0 0 0], ...
    'Marker', 'o',...
    'LineStyle','none', ...
    'MarkerSize', 10, ...
    'DisplayName','Active');

plot(x(headList(find(headList>0))),y(headList(find(headList>0))), ... 
    'MarkerFaceColor',[0.2 1 0.2], ...
    'MarkerEdgeColor',[0.3 0.3 0.3], ...
    'Marker','s',...
    'LineStyle','none', ...
    'MarkerSize',15, ...
    'DisplayName','Head');

%plot node topology
plot(0,0,'^', ... 
    'MarkerSize',20, ...
    'MarkerFaceColor','r',...
    'MarkerEdgeColor','k',...
    'DisplayName','BS',...
    'LineStyle','none'); %Base Station

%format
%axis([-(radius/1.7+80) (radius/1.7+80) -(radius/1.7+150) (radius/1.7+80)]);
%axis([-110 90 -50 150]);
axis([-(radius/1.2) (radius/1.2) -(radius/1.2) (radius/1.2)]);
%Write Interpretation

str = sprintf('Machines=%d, \\eta=0.48',totalNodes);

%title('Machines=50, \eta=0.477');
%plot X-Y Axis
temp = -(radius+10):0.01:(radius+10);
%plot (temp,zeros(1,length(temp)),'k');     
%plot (zeros(1,length(temp)),temp),'k';
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'YColor','w');
set(gca,'XColor','w');
set(gca,'XTickLabel',['']);
set(gca,'YTickLabel',['']);
% xlabel({'x-axis(m)'});
% ylabel('y-axis(m)');
hold off;
%legend1 = legend('show','Orientation','horizontal');
%set(legend1,...
%    'Position',[0.23590909090909 0.854000000000001 0.564545454545455 0.0025]);
end
