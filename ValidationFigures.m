clear all
close all

fileloc = 'U:\Research\Projects\ihbi\btm\btm_general\BTM Folder Structure\Mark_Allenby\05_Research\05_Research\Aneurysm_Data_DevelopmentCode\IXI_120exp_20052020';
% IXI_178_24052020  IXI_60exp_24052020  IXI_90exp_24052020  IXI_120exp_20052020
load(fileloc)
REP1.UIAcounter = UIAcounter;
REP1.data_subpl = data_subpl;
REP1.fitresult = fitresult;
REP1.time2run = time2run;
REP1.outlierProps = outlierProps;
REP1.KwPredictor = KwPredictor;
REP1.PathLength = PathLength;
REP1.MaskVol = MaskVol;

fileloc = 'U:\Research\Projects\ihbi\btm\btm_general\BTM Folder Structure\Mark_Allenby\05_Research\Aneurysm_Data_DevelopmentCode\MIDAS_120exp_20052020';
% MIDAS_178_24052020  MIDAS_60exp_20052020 MIDAS_90exp_20052020  MIDAS_120exp_20052020
load(fileloc)
REP2.UIAcounter = UIAcounter;
REP2.data_subpl = data_subpl;
REP2.fitresult = fitresult;
REP2.time2run = time2run;
REP2.outlierProps = outlierProps;
REP2.KwPredictor = KwPredictor;
REP2.PathLength = PathLength;
REP2.MaskVol = MaskVol;

fileloc = 'U:\Research\Projects\ihbi\btm\btm_general\BTM Folder Structure\Mark_Allenby\05_Research\Aneurysm_Data_DevelopmentCode\RBWH_120exp_25052020.xlsx';
% RBWH_178_15052020  RBWH_60exp_27052020  RBWH_90exp_21052020  RBWH_120exp_25052020
RBWHtab = readtable(fileloc);
RBWH.Names = join([table2array(RBWHtab(:,1)), cellstr(num2str(table2array(RBWHtab(:,2))))]);
RBWH.UIAcounter = table2array(RBWHtab(:,7));
RBWH.time2run = (table2array(RBWHtab(:,11)));
RBWH.KwPredictor = (table2array(RBWHtab(:,12)));
RBWH.PathLength = (table2array(RBWHtab(:,14)));
RBWH.MaskVol = (table2array(RBWHtab(:,13)));
RBWH.fitresult = (table2array(RBWHtab(:,16:24)));

fileloc = 'U:\Research\Projects\ihbi\btm\btm_general\BTM Folder Structure\Mark_Allenby\05_Research\Aneurysm_Data_DevelopmentCode\RBWH_120exp_Healthy_25052020.xlsx';
% RBWH_178_Healthy_15052020  RBWH_60exp_Healthy_27052020  RBWH_90exp_Healthy_15052020  RBWH_120exp_Healthy_25052020
RBWHhtab = readtable(fileloc);
RBWHh.Names = join([table2array(RBWHhtab(:,1)), cellstr(num2str(table2array(RBWHhtab(:,2))))]);
RBWHh.UIAcounter = table2array(RBWHhtab(:,7));
RBWHh.time2run = (table2array(RBWHhtab(:,11)));
RBWHh.KwPredictor = (table2array(RBWHhtab(:,12)));
RBWHh.PathLength = (table2array(RBWHhtab(:,14)));
RBWHh.MaskVol = (table2array(RBWHhtab(:,13)));
RBWHh.fitresult = (table2array(RBWHhtab(:,16:24)));

RBWHttab = [RBWHtab;RBWHhtab];
RBWHt.UIAcounter = table2array(RBWHttab(:,7));
RBWHt.correct = [RBWH.UIAcounter==1; RBWHh.UIAcounter==0];
RBWHt.correct(isnan(RBWHt.UIAcounter)) = [];
RBWHttab(isnan(RBWHt.UIAcounter),:) = [];
RBWHt.Names = join([table2array(RBWHttab(:,1)), cellstr(num2str(table2array(RBWHttab(:,2))))]);
RBWHt.UIAcounter = table2array(RBWHttab(:,7));
RBWHt.time2run = (table2array(RBWHttab(:,11)));
RBWHt.KwPredictor = (table2array(RBWHttab(:,12)));
RBWHt.PathLength = (table2array(RBWHttab(:,14)));
RBWHt.MaskVol = (table2array(RBWHttab(:,13)));
RBWHt.fitresult = (table2array(RBWHttab(:,16:24)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scatterplot Mask Length and Mask Volume and UIA Detection

figure
scatter(REP1.PathLength(REP1.UIAcounter>0),REP1.MaskVol(REP1.UIAcounter>0),'ko','MarkerEdgeAlpha', 1,'LineWidth',1)
hold on
scatter(RBWH.PathLength(RBWH.UIAcounter==0),RBWH.MaskVol(RBWH.UIAcounter==0),'ro','MarkerEdgeAlpha', 1,'LineWidth',1)
scatter(REP1.PathLength,REP1.MaskVol,'bo','filled','MarkerFaceAlpha', 0.25)
scatter(REP2.PathLength,REP2.MaskVol,'ro','filled','MarkerFaceAlpha', 0.25)
scatter(RBWHt.PathLength,RBWHt.MaskVol,'co','filled','MarkerFaceAlpha', 1,'MarkerEdgeColor',[0.2 0.4 0.7],'LineWidth',0.5)
scatter(RBWH.PathLength(RBWH.UIAcounter==0),RBWH.MaskVol(RBWH.UIAcounter==0),'co','filled','MarkerFaceAlpha', 1,'MarkerEdgeColor','r','LineWidth',1)
scatter(REP1.PathLength(REP1.UIAcounter>0),REP1.MaskVol(REP1.UIAcounter>0),'ko','MarkerEdgeAlpha', 1)
scatter(RBWHh.PathLength(RBWHh.UIAcounter>0),RBWHh.MaskVol(RBWHh.UIAcounter>0),'co','filled','MarkerFaceAlpha', 1,'MarkerEdgeColor','k','LineWidth',1)
scatter(REP2.PathLength(REP2.UIAcounter>0),REP2.MaskVol(REP2.UIAcounter>0),'ko','MarkerEdgeAlpha', 1)
scatter(RBWHt.PathLength,RBWHt.MaskVol,'co','filled','MarkerFaceAlpha', 1,'MarkerEdgeColor',[0.2 0.4 0.7],'LineWidth',0.5)
scatter(RBWH.PathLength(RBWH.UIAcounter==0),RBWH.MaskVol(RBWH.UIAcounter==0),'co','filled','MarkerFaceAlpha', 1,'MarkerEdgeColor','r','LineWidth',1)
scatter(RBWHh.PathLength(RBWHh.UIAcounter>0),RBWHh.MaskVol(RBWHh.UIAcounter>0),'co','filled','MarkerFaceAlpha', 1,'MarkerEdgeColor','k','LineWidth',1)
ylabel('Mask Volume (mm^3)')
xlabel('Mask Length (mm)')
leg = legend([strcat('Total specificity:',{' '},num2str(round(sum([REP1.UIAcounter==0,REP2.UIAcounter==0,RBWHh.UIAcounter'==0])/length([REP1.UIAcounter,REP2.UIAcounter,RBWHh.UIAcounter'])*100)),'%',...
    {' '},'(n=',num2str(size(REP1.UIAcounter,2)+size(REP2.UIAcounter,2)+size(RBWH.UIAcounter,1)),')'),...
    strcat('RBWH sensitivity:',{' '},num2str(round(sum(RBWH.UIAcounter==1)/sum(~isnan(RBWH.UIAcounter))*100)),'%',{' '},'(n=',num2str(size(RBWHh.UIAcounter,1)),')'),...
    strcat('IXI:             n=',num2str(size(REP1.UIAcounter,2)),', spec:',{' '},num2str(round((sum(REP1.UIAcounter==0)/size(REP1.UIAcounter,2))*100)),'%'), ...
    strcat('MIDAS:      n=',num2str(size(REP2.UIAcounter,2)),', spec:',{' '},num2str(round((sum(REP2.UIAcounter==0)/size(REP2.UIAcounter,2))*100)),'%'),...
    strcat('RBWH:      n=',num2str(length(RBWHt.UIAcounter)),', spec:',{' '},num2str(round(((sum(RBWHh.UIAcounter==0)/sum(~isnan(RBWHh.UIAcounter))))*100)),'%')],'Location','northwest','Orientation','Vertical');
title(leg,'TOF MRA Validation')
set(gca,'xscale','log','yscale','log','XLim',[3*10^1 6*10^2],'YLim',[3*10^2 10^5])
set(gca,'XTick',[30 100 600],'YTick',[3*10^2 get(gca,'YTick')])
set(gca,'YTickLabel',[{'300'} {'10^3'} {'10^4'} {'10^5'}]) 
set(gcf, 'Position',  [100, 100, 500, 475])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Histogram Mask Length

clear Count1 Count2 CountR
res = 0.1;
bins = 500; %20
binedges = logspace(log10(30),log10(600),bins+1);
bincentres = logspace(log10(30),log10(600),(bins+1)*2-1);
bincentres = bincentres(2:2:(length(bincentres)-1));
for i=1:(length(binedges)-1)
   Count1(i)  = sum((REP1.PathLength>binedges(i)).*(REP1.PathLength<binedges(i+1)));
   Count2(i)  = sum((REP2.PathLength>binedges(i)).*(REP2.PathLength<binedges(i+1)));
   CountR(i)  = sum((RBWHt.PathLength>binedges(i)).*(RBWHt.PathLength<binedges(i+1)))';
end
Smooth1 = smooth(bincentres,Count1,res,'lowess');
Smooth2 = smooth(bincentres,Count2,res,'lowess');
SmoothR = smooth(bincentres,CountR,res,'lowess');
figure
hold on
scatter(bincentres,Smooth1/sum(Smooth1),10,'k','MarkerFaceColor','b','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.1)
scatter(bincentres,Smooth2/sum(Smooth2),10,'k','MarkerFaceColor','r','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.1)
scatter(bincentres,SmoothR/sum(SmoothR),10,'k','MarkerFaceColor','c','MarkerFaceAlpha',1,'MarkerEdgeAlpha',0.1)
set(gca,'xscale','log','XLim',[3*10^1 6*10^2],'XTick',[30 100 600],'YTick',[])
set(gcf, 'Position',  [100, 100, 500, 100])
xlabel('Mask Length (mm)')
leg = legend(strcat(num2str(round(mean(REP1.PathLength))),'\pm',num2str(round(std(REP1.PathLength))),'mm'),...
    strcat(num2str(round(mean(REP2.PathLength))),'\pm',num2str(round(std(REP2.PathLength))),'mm'),...
    strcat(num2str(round(mean(RBWHt.PathLength))),'\pm',num2str(round(std(RBWHt.PathLength))),'mm'),...
    'Location','northwest','Orientation','Vertical');

% figure
% hold on
% histogram('BinEdges',binedges,'BinCounts',Count1,'Normalization','probability','FaceColor','b','FaceAlpha',0.25)
% histogram('BinEdges',binedges,'BinCounts',Count2,'Normalization','probability','FaceColor','r','FaceAlpha',0.25)
% histogram('BinEdges',binedges,'BinCounts',CountR,'Normalization','probability','FaceColor','c','FaceAlpha',0.5)
% set(gca,'xscale','log','XLim',[3*10^1 6*10^2],'XTick',[30 100 600],'YTick',[])
% set(gcf, 'Position',  [100, 100, 500, 100])
% xlabel('Mask Length (mm)')
% leg = legend(strcat(num2str(round(mean(REP1.PathLength))),'\pm',num2str(round(std(REP1.PathLength))),'mm'),...
%     strcat(num2str(round(mean(REP2.PathLength))),'\pm',num2str(round(std(REP2.PathLength))),'mm'),...
%     strcat(num2str(round(mean(RBWHt.PathLength))),'\pm',num2str(round(std(RBWHt.PathLength))),'mm'),...
%     'Location','northwest','Orientation','Vertical');

% pd1 = fitdist(log10(REP1.PathLength)','Normal');
% y1 = pdf(pd1,log10(bincentres));
% pd2 = fitdist(log10(REP2.PathLength)','Normal');
% y2 = pdf(pd2,log10(bincentres));
% pdR = fitdist(log10(RBWHt.PathLength),'Normal');
% yR = pdf(pdR,log10(bincentres));
% figure
% hold on
% area(bincentres,y1/sum(y1),'EdgeColor','b','FaceColor','b','FaceAlpha',0.25)
% area(bincentres,y2/sum(y2),'EdgeColor','r','FaceColor','r','FaceAlpha',0.25)
% area(bincentres,yR/sum(yR),'EdgeColor','k','FaceColor','c','FaceAlpha',1)
% set(gca,'xscale','log','XLim',[3*10^1 6*10^2],'XTick',[30 100 600],'YTick',[])
% set(gcf, 'Position',  [100, 100, 500, 100])
% xlabel('Mask Length (mm)')
% leg = legend(strcat(num2str(round(mean(REP1.PathLength))),'\pm',num2str(round(std(REP1.PathLength))),'mm'),...
%     strcat(num2str(round(mean(REP2.PathLength))),'\pm',num2str(round(std(REP2.PathLength))),'mm'),...
%     strcat(num2str(round(mean(RBWHt.PathLength))),'\pm',num2str(round(std(RBWHt.PathLength))),'mm'),...
%     'Location','northwest','Orientation','Vertical');

[h,p] = ttest2(RBWHt.PathLength,REP1.PathLength,'Vartype','equal')
[h,p] = ttest2(RBWHt.PathLength,REP2.PathLength,'Vartype','equal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Histogram Mask Volume

clear Count1 Count2 CountR
binedges = logspace(log10(300),log10(10^5),bins+1);
bincentres = logspace(log10(300),log10(10^5),(bins+1)*2-1);
bincentres = bincentres(2:2:(length(bincentres)-1));
for i=1:(length(binedges)-1)
   Count1(i)  = sum((REP1.MaskVol>binedges(i)).*(REP1.MaskVol<binedges(i+1)));
   Count2(i)  = sum((REP2.MaskVol>binedges(i)).*(REP2.MaskVol<binedges(i+1)));
   CountR(i)  = sum((RBWHt.MaskVol>binedges(i)).*(RBWHt.MaskVol<binedges(i+1)))';
end
Smooth1 = smooth(bincentres,Count1,res,'lowess');
Smooth2 = smooth(bincentres,Count2,res,'lowess');
SmoothR = smooth(bincentres,CountR,res,'lowess');
figure
hold on
scatter(bincentres,Smooth1/sum(Smooth1),10,'k','MarkerFaceColor','b','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.1)
scatter(bincentres,Smooth2/sum(Smooth2),10,'k','MarkerFaceColor','r','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.1)
scatter(bincentres,SmoothR/sum(SmoothR),10,'k','MarkerFaceColor','c','MarkerFaceAlpha',1,'MarkerEdgeAlpha',0.1)
 set(gca,'xscale','log','XLim',[3*10^2 10^5])
 set(gca,'XTick',[3*10^2 get(gca,'XTick')],'XTickLabel',[{'300'} {'10^3'} {'10^4'} {'10^5'}],'YTick',[]) 
 set(gcf, 'Position',  [100, 100, 500, 100])
 xlabel('Mask Volume (mm^3)')
leg = legend(strcat(num2str(round(mean(REP1.MaskVol))),'\pm',num2str(round(std(REP1.MaskVol))),'mm^3'),...
    strcat(num2str(round(mean(REP2.MaskVol))),'\pm',num2str(round(std(REP2.MaskVol))),'mm^3'),...
    strcat(num2str(round(mean(RBWHt.MaskVol))),'\pm',num2str(round(std(RBWHt.MaskVol))),'mm^3'),...
    'Location','northeast','Orientation','Vertical');

% pd1 = fitdist(log10(REP1.MaskVol)','Normal');
% y1 = pdf(pd1,log10(bincentres));
% pd2 = fitdist(log10(REP2.MaskVol)','Normal');
% y2 = pdf(pd2,log10(bincentres));
% pdR = fitdist(log10(RBWHt.MaskVol),'Normal');
% yR = pdf(pdR,log10(bincentres));
% figure
% hold on
% % area(bincentres,y1/sum(y1),'EdgeColor','b','FaceColor','b','FaceAlpha',0.25)
% % area(bincentres,y2/sum(y2),'EdgeColor','r','FaceColor','r','FaceAlpha',0.25)
% % area(bincentres,yR/sum(yR),'EdgeColor','k','FaceColor','c','FaceAlpha',1)
% scatter(bincentres,y1/sum(y1),10,'k','MarkerFaceColor','b','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.1)
% scatter(bincentres,y2/sum(y2),10,'k','MarkerFaceColor','r','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.1)
% scatter(bincentres,yR/sum(yR),10,'k','MarkerFaceColor','c','MarkerFaceAlpha',1,'MarkerEdgeAlpha',0.1)
% set(gca,'xscale','log','XLim',[3*10^2 10^5])
% set(gca,'XTick',[3*10^2 get(gca,'XTick')],'XTickLabel',[{'300'} {'10^3'} {'10^4'} {'10^5'}],'YTick',[]) 
% set(gcf, 'Position',  [100, 100, 500, 100])
% xlabel('Mask Volume (mm^3)')
% leg = legend(strcat(num2str(round(mean(REP1.MaskVol))),'\pm',num2str(round(std(REP1.MaskVol))),'mm^3'),...
%     strcat(num2str(round(mean(REP2.MaskVol))),'\pm',num2str(round(std(REP2.MaskVol))),'mm^3'),...
%     strcat(num2str(round(mean(RBWHt.MaskVol))),'\pm',num2str(round(std(RBWHt.MaskVol))),'mm^3'),...
%     'Location','northeast','Orientation','Vertical');

[h,p] = ttest2(RBWHt.MaskVol,REP1.MaskVol,'Vartype','equal')
[h,p] = ttest2(RBWHt.MaskVol,REP2.MaskVol,'Vartype','equal')

% [~,edges] = histcounts(log10(REP1.PathLength));
% h = histogram(REP1.PathLength,10.^edges,'Normalization','pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Median Mask Regression

REP1.fitresult(REP1.UIAcounter>0,:) = [];
REP2.fitresult(REP2.UIAcounter>0,:) = [];
RBWHt.fitresult(RBWHt.correct>0,:) = [];
x = repmat(linspace(0,200,500),500,1);
y = repmat(linspace(0.5,2.5,500)',1,500);

figure
hold on
xlim([0 200]) %325
ylim([0 3])
zlim([0 2]) %8
view(135,30);
for i = 1:size(REP1.fitresult,1)
    Pval1 = REP1.fitresult(i,:);
    REP1z(:,:,i) = Pval1(1)+Pval1(2)*x+Pval1(3)*y+Pval1(4)*x.^2+Pval1(5)*x.*y+Pval1(6)*x.^3+...
        Pval1(7)*x.^2.*y+Pval1(8)*x.^4+Pval1(9)*x.^3.*y;
end
REP1zm = median(REP1z,3);
for i = 1:size(REP2.fitresult,1)
    Pval1 = REP2.fitresult(i,:);
    REP2z(:,:,i) = Pval1(1)+Pval1(2)*x+Pval1(3)*y+Pval1(4)*x.^2+Pval1(5)*x.*y+Pval1(6)*x.^3+...
        Pval1(7)*x.^2.*y+Pval1(8)*x.^4+Pval1(9)*x.^3.*y;
end
REP2zm = median(REP2z,3);
for i = 1:size(RBWHt.fitresult,1)
    Pval1 = RBWHt.fitresult(i,:);
    RBWHz(:,:,i) = Pval1(1)+Pval1(2)*x+Pval1(3)*y+Pval1(4)*x.^2+Pval1(5)*x.*y+Pval1(6)*x.^3+...
        Pval1(7)*x.^2.*y+Pval1(8)*x.^4+Pval1(9)*x.^3.*y;
end
RBWHzm = median(RBWHz,3);
surf(x,y,REP1zm,'FaceColor',[0 0 0.5],'FaceAlpha',0.1,'EdgeColor','none')
surf(x,y,REP2zm,'FaceColor',[0.5 0 0],'FaceAlpha',0.1,'EdgeColor','none')
surf(x,y,RBWHzm,'FaceColor',[0.25 1 1],'FaceAlpha',0.2,'EdgeColor','none')
plot3(x(1,:),y(1,:),REP1zm(1,:),'b-')
plot3(x(:,1),y(:,1),REP1zm(:,1),'b-')
for i = 1:size(x,1)
    [zint(i),yint1(i)] = min(abs(REP1zm(:,i)));
end
ind1 = (zint<0.02).*(1:size(zint,2));
ind1(ind1==0) = [];
plot3(x(1,ind1),y(yint1(ind1),1),zeros(1,size(ind1,2)),'b-')
plot3(x(1,:),y(1,:),REP2zm(1,:),'r-')
plot3(x(:,1),y(:,1),REP2zm(:,1),'r-')
for i = 1:size(x,1)
    [zint(i),yint2(i)] = min(abs(REP2zm(:,i)));
end
ind2 = (zint<0.02).*(1:size(zint,2));
ind2(ind2==0) = [];
plot3(x(1,ind2),y(yint2(ind2),1),zeros(1,size(ind2,2)),'r-')
plot3(x(1,:),y(1,:),RBWHzm(1,:),'Color',[0 0.5 0.5])
plot3(x(:,1),y(:,1),RBWHzm(:,1),'Color',[0 0.5 0.5])
for i = 1:size(x,1)
    [zint(i),yint3(i)] = min(abs(RBWHzm(:,i)));
end
ind3 = (zint<0.02).*(1:size(zint,2));
ind3(ind3==0) = [];
plot3(x(1,ind3),y(yint3(ind3),1),zeros(1,size(ind3,2)),'Color',[0 0.5 0.5])
xlabel('Base Distance {\it d_b} (mm)')
ylabel('Edge Distance {\it d_e} (mm)')
zlabel('Centreline Distance {\it d_c} (mm)')

REP1zOm = REP1zm;
REP1zOm(REP1zOm<0)=0;
REP2zOm = REP2zm;
REP2zOm(REP2zOm<0)=0;
RBWHzOm = RBWHzm;
RBWHzOm(RBWHzOm<0)=0;
R1pDiff = (abs(REP1zOm - RBWHzOm));
R1pDiff(R1pDiff==0) = [];
R2pDiff = (abs(REP2zOm - RBWHzOm));
R2pDiff(R2pDiff==0) = [];
 leg = legend([strcat('IXI',{'          '},char(949),{' '},{'='},{' '},num2str(round(mean(R1pDiff)*100)/100),'\pm',num2str(round(std(R1pDiff)*100)/100),'mm'), ...
     strcat('MIDAS',{'   '},char(949),{' '},{'='},{' '},num2str(round(mean(R2pDiff)*100)/100),'\pm',num2str(round(std(R2pDiff)*100)/100),'mm'),...
     strcat('RBWH',{' '})],'Location','northwest','Orientation','Vertical');
title(leg,'Median Voxel Regression')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list = (REP2.UIAcounter>0).*(1:length(REP2.UIAcounter));
list(list==0) = [];
for i = 1:length(list)
    clear DensityTemp
    for j = 1:length(REP2.data_subpl{list(i)})
        PixelCoord{j} = REP2.data_subpl{list(i)}(j).PixelList;
     DensityTemp(j) = size(PixelCoord{j},1)/prod(max(PixelCoord{j}) - min(PixelCoord{j}));
    end
    if j>1
    [~, b] = max(DensityTemp);
    PixSize(i) = size(PixelCoord{b},1);
    MinSlice(i) = min(PixelCoord{b}(:,3));
    else
        PixSize(i) = size(PixelCoord{1},1);
        MinSlice(i) = min(PixelCoord{1}(:,3));
    end
    Density(i,:) = [list(i) max(DensityTemp)];
end
[Dense b] = sort(Density(:,2),'descend');
LostImages = unique([PoorImages DeleteList]);
Ordered{1} = list(b)';
Ordered{2} = Ordered{1} + sum(Ordered{1}>=LostImages,2);
i=2;
while sum(sum(Ordered{i}>=LostImages,2)-sum(Ordered{i-1}>=LostImages,2))>0
i=i+1;
Ordered{i} = Ordered{i-1} + (sum(Ordered{i-1}>=LostImages,2)-sum(Ordered{i-2}>=LostImages,2));
end
[Ordered{end} Dense PixSize(b)' MinSlice(b)']
size(LostImages)


list = (REP1.UIAcounter>0).*(1:length(REP1.UIAcounter));
list(list==0) = [];
for i = 1:length(list)
    clear DensityTemp
    for j = 1:length(REP1.data_subpl{list(i)})
        PixelCoord{j} = REP1.data_subpl{list(i)}(j).PixelList;
        DensityTemp(j) = size(PixelCoord{j},1)/prod(max(PixelCoord{j}) - min(PixelCoord{j}));
    end
    if j>1
    [~, b] = max(DensityTemp);
    PixSize(i) = size(PixelCoord{b},1);
    MinSlice(i) = min(PixelCoord{b}(:,3));
    else
        PixSize(i) = size(PixelCoord{1},1);
        MinSlice(i) = min(PixelCoord{1}(:,3));
    end
    Density(i,:) = [list(i) max(DensityTemp)];
end
[Dense b] = sort(Density(:,2),'descend');
LostImages = unique([PoorImages DeleteList]);
Ordered{1} = list(b)';
Ordered{2} = Ordered{1} + sum(Ordered{1}>=LostImages,2);
i=2;
while sum(sum(Ordered{i}>=LostImages,2)-sum(Ordered{i-1}>=LostImages,2))>0
i=i+1;
Ordered{i} = Ordered{i-1} + (sum(Ordered{i-1}>=LostImages,2)-sum(Ordered{i-2}>=LostImages,2));
end
[Ordered{end} Dense PixSize(b)' MinSlice(b)']
size(LostImages)
