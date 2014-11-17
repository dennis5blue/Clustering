objGo = [0.7310,0.7315,inf,inf,0.7429,inf,0.8498,inf,0.8378,inf,0.9447,inf, ...
         0.7311,0.7430,0.8379,0.7314,inf,inf,inf,0.8015,inf,inf,0.8964,inf, ...
         0.7347,0.8415,0.7350,0.7352,inf,0.8561,inf,inf,0.7937,inf,0.9,0.9485,1.2918,1.1370, ...
         0.7550,0.7669,0.7553,0.7555,inf,0.8622,inf,inf,0.8140,inf,0.8253,0.8738,1.2172,1.0624];
     
objUpdate = [inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,inf, ...
             inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,inf, ...
             inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,1.2918,1.1370, ...
             1.1370,1.1370,1.1370,1.1370,1.1370,1.1370,1.1370,1.1370,1.1370,1.1370,1.1370,1.1370,1.1370,1.0624];
iter = [1:52];

str_objGo = sprintf('Lower bound');
plot(iter,objGo,'x','LineWidth',2,'Color','b','DisplayName',str_objGo); hold on;
str_objUpdate = sprintf('Upper bound');
plot(iter,objUpdate,'-o','LineWidth',2,'Color','r','DisplayName',str_objUpdate); hold on;

axis([-inf,inf,0.7,1.3]);
xlabel('Iteration (run)','FontSize',13);
ylabel('Tier-2 Power (Watt)','FontSize',13);
legend('show','location','NorthWest');
grid on;