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

path = {'./data/CS_test_DMCP.out'};

strtitle2={'SA'};
mapString = ['./data/20cam_r500_map.out'];
vecDirection = dlmread('./data/20cam_r500_dir.out');
%path=cellstr(dpath);
%strtitle2=cellstr(strtitle);
for ii=1:length(path)
% 
hFig = figure(ii);
set(hFig, 'Position', [1 1 550 500 ]);
%subplot(2,2,ii); 
%xp=0.07+0.31*(ii-1);

%subplot('Position',[ xp 0.1 0.24 0.8]);
%Show up later: Gij = dBToRaw(-(131.1 + 42.81*log10( (D2.^0.5)/1000 )));
[totalNodes radius x y Gij Gib]= parse_Map_v1(mapString);
powerBound = 1e-6;
[ maxChNum powerMax C2WT payoffs SAFac TimeMs1st TimeMs2nd EnergyJ1st EnergyJ2nd clusterStru headList TxPower  ] = parse_Stru_v2 ...
    (path{ii},totalNodes);
headName=headList;

radius = 500;
densityScale = 1000;%km square
density = totalNodes / (pi*(radius/densityScale)^2) ;
%Aboves are basic parameter input

hold on;
%circle(-10,55,radius/1.7);
circle(0,0,radius/1.2);

%plot (x,y,'o','MarkerSize',4);%all nodes

B = find(sum(clusterStru));
supSet = setdiff(B, headList);


Binary_Unsupset=ones(1,totalNodes)-sum(clusterStru);
index_unsupset=find(Binary_Unsupset);

%circle(-10,55,radius/1.7);
circle(0,0,radius/1.2);
%grid on;
clusterSize=zeros(1,maxChNum);

  firstTierX(1)=0;
  firstTiery(1)=0;
%connect the head and serves nodes
for i=1:maxChNum
  if(headName(i)==0)continue;end
  firstTierX(2)=x(headName(i));
  firstTiery(2)=y(headName(i));
  hLine=plot(firstTierX,firstTiery,'k-','LineWidth',1.7);
  set(get(get(hLine,'Annotation'),'LegendInformation'),...
  'IconDisplayStyle','off');
  for j=1:totalNodes
  if clusterStru(i,j)==1
     X(1) = x(headName(i));
     Y(1) = y(headName(i));
     X(2) = x(j);
     Y(2) = y(j);
     theta = vecDirection(j)+pi/2;
     r = 20;
     u = x(j) + r * cos(theta); % convert polar (theta,r) to cartesian
     v = y(j) + r * sin(theta);
     
     arrow([x(j) y(j)],[u v],'BaseAngle',30,'length',15,'width',0.4);
     if i==1
         hLine = plot(X,Y,'k:','Color','r','LineWidth',3);
         set(get(get(hLine,'Annotation'),'LegendInformation'),...
         'IconDisplayStyle','off');
     elseif i==2
         hLine = plot(X,Y,'k:','Color','g','LineWidth',3);
         set(get(get(hLine,'Annotation'),'LegendInformation'),...
         'IconDisplayStyle','off');
     elseif i==3
         hLine = plot(X,Y,'k:','Color','b','LineWidth',3);
         set(get(get(hLine,'Annotation'),'LegendInformation'),...
         'IconDisplayStyle','off');
     elseif i==4
         hLine = plot(X,Y,'k:','Color','c','LineWidth',3);
         set(get(get(hLine,'Annotation'),'LegendInformation'),...
         'IconDisplayStyle','off');
     end
     clusterSize(i) = clusterSize(i)+1;
     
     %str = sprintf('%d',j);
     %text(x(j)-10,y(j)+10,str,'FontSize',10);
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

%%
%Compute the SINR below

 [Rate SINR Inter]  = calSINR_4All( maxChNum, totalNodes, clusterStru,  clusterSize,  headName, TxPower,  Gij  );
%Add text
for i=1:maxChNum
  
  if(headName(i)==0)continue;end
   for j=1:totalNodes
    if (headName(i)==j && clusterStru(i,j)==1)
     %str = sprintf('H %d', headName(i));
     %text(x(j)+10,y(j)-10,str,'FontSize',12);
     %str = sprintf('C%d',i);
     %text(x(j)+10,y(j)-10,str,'FontSize',10);
  
    elseif clusterStru(i,j)==1
     %str = sprintf('%.5f',SINR(j));
     %text(x(j),y(j),str,'FontSize',12);
    end
  end
end
aa = plot(x(index_unsupset),y(index_unsupset),'LineWidth',3,'MarkerFaceColor',[0.2 1 0],'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',10,...
    'Marker','X',...
    'LineStyle','none', 'DisplayName','Unselected Machine');
bb = plot(x(supSet),y(supSet),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'Marker','o',...
   'LineStyle','none', 'MarkerSize', 5, 'DisplayName','Member');
cc = plot(x(headList(find(headList>0))),y(headList(find(headList>0))),'MarkerFaceColor',[0.3 1 0.3],'MarkerEdgeColor',[0.3 0.3 0.3],'Marker','s',...
   'LineStyle','none','MarkerSize',7, 'DisplayName','Head');

%plot node topology
dd = plot(0,0,'^','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','BS','LineStyle','none');%Base Station





%format
%axis([-(radius+10) (radius+10) -(radius+10) (radius+10)]);
%axis([-110 90 -50 150]);
axis([-(radius/1.2) (radius/1.2) -(radius/1.2) (radius/1.2)]);
%Write Interpretation

%str = sprintf('UL Structure(N%d-HN%d-Served:%d-RT%.3f(bps/Hz)-Density(%.3f/km^2)-PowerMax %.3f(W/Hz). SAloop %d in %.2f seconds)',totalNodes, maxChNum,payoffs, C2W, density, powerBound,SAFac,computingTime);
strtitle2{ii};
str = sprintf('Machines=50, \lambda=0.8 \eta=0.477:%s ',strtitle2{ii});

title('Cameras=20');
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
legend('show','Orientation','horizontal');legend([aa bb cc dd]);
legend1 = legend('show');
set(legend1,...
    'Position',[0.65590909090909 0.754000000000001 0.194545454545455 0.1025]);



end

