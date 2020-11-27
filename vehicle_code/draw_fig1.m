
clear;clc;close all;
load data1;
c=real(cmds);w=real(wmds);g=real(gps);

Jan_O=g(:,1);Feb_O=g(:,2);Mar_O=g(:,3);Apr_O=g(:,4);May_O=g(:,5);
%Jun_O=g(:,6);Jul_O=g(:,7);Aug_0=g(:,8);
Temp_O = [Jan_O, Feb_O, Mar_O, Apr_O, May_O];
position_O = 1.2:1:5.2; 
% Define position for 12 Month_O boxplots 
box_O = boxplot(Temp_O,'colors','r','positions',position_O,'width',0.18,'sym',' ','whisker',2);hold on;
set(gca,'XTickLabel',{' '})  % Erase xlabels  
hold on  % Keep the Month_O boxplots on figure overlap the Month_S boxplots  

Jan_Z=w(:,1);Feb_Z=w(:,2);Mar_Z=w(:,3);Apr_Z=w(:,4);May_Z=w(:,5);
%Jun_Z=w(:,6);Jul_Z=w(:,7);Aug_Z=w(:,8);
Temp_Z = [Jan_Z, Feb_Z, Mar_Z, Apr_Z, May_Z];
position_Z = 1.6:1:5.6;  % Define position for 12 Month_S boxplots 
box_Z = boxplot(Temp_Z,'colors','b','positions',position_Z,'width',0.18,'sym',' ','whisker',2);  

Jan_S=c(:,1);Feb_S=c(:,2);Mar_S=c(:,3);Apr_S=c(:,4);May_S=c(:,5);
%Jun_S=c(:,6);Jul_S=c(:,7);Aug_S=c(:,8);
Temp_S = [Jan_S, Feb_S, Mar_S, Apr_S, May_S];
position_S = 1.4:1:5.4;  % Define position for 12 Month_S boxplots 
box_S = boxplot(Temp_S,'colors','g','positions',position_S,'width',0.18,'sym',' ','whisker',2);

% legend('MDS-MAP','MDS-MAP(P)','SMDS','SMDS(P)','GPS');
box_vars=findall(gca,'Tag','Box');
hLegend = legend(box_vars([15,1,10]), {'GPS','Classical MDS','Optimized MDS'});

xlabel('$\alpha$','interpreter','latex','fontsize',16);% M
ylabel('Error (m)');
set(gca,'XLim',[0.8 6]);
set(gca,'YLim',[0 0.6]);

