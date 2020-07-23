function [data_subpl,fitresult1,blanks,blankp,blankd,blankdiam,outlierProps] = AnomalyDetection(Images,data_pl,blanks,blanki,blankdiam,blankh,...
    blankd,blankp,Matched,ConInt,minClustRate,SmallLinks,INFO,BottomBranch,RefineIA,ConIntRefine,Dclose,RRadMult)

if RefineIA == 1 %if refinement occurs, needs two separate outlier analyses. 
    ConInt2 = ConInt;
    ConInt = ConIntRefine;
end

data_pl = data_pl(logical((data_pl(:,3)>size(Images,3)*0.1).*(data_pl(:,3)<size(Images,3)*0.95)),:);
[fitresult1] = fit([blankp(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3))),...
    blankd(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3)))],...
    blankdiam(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3))),'poly41');
[pci py] = predint(fitresult1,[blankp(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3))) ...
    blankd(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3)))],ConInt,'observation','off');
ry = blankdiam(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3)));
SubSec = data_pl(logical(ry>pci(:,2)),:);

blankSub = zeros(size(Images));
blankSub(sub2ind(size(Images),SubSec(:,2), SubSec(:,1), SubSec(:,3))) = 1;
data_subc = bwconncomp(blankSub,6);
data_subpl = regionprops(data_subc, 'PixelList');
Voxels_sub = cell2mat(cellfun(@length,cellfun(@transpose,{data_subpl.PixelList},'UniformOutput',false),'uni',false));

try
    i = 1;
    while i <= size(data_subpl,1)
        OutlierProps = mean(fitresult1(blankp(sub2ind(size(Images),data_subpl(i).PixelList(:,2),data_subpl(i).PixelList(:,1),data_subpl(i).PixelList(:,3))),...
            blankd(sub2ind(size(Images),data_subpl(i).PixelList(:,2),data_subpl(i).PixelList(:,1),data_subpl(i).PixelList(:,3)))));
        OutlierProps(OutlierProps<0) = 0;
        if Voxels_sub(i)<(OutlierProps*minClustRate(1)+minClustRate(2))
            data_subpl(i) = [];
            Voxels_sub(i) = [];
        else
            i = i+1;
        end
    end
end

if RefineIA == 1 %the start of opional mask refinement/culling
    SubListInit = vertcat(data_subpl.PixelList);
    for j = 1:length(data_subpl)
        hold on
        scatter3(mean(data_subpl(j).PixelList(:,1))*INFO.PixelSpacing(1),mean(data_subpl(j).PixelList(:,2))*INFO.PixelSpacing(2),mean(data_subpl(j).PixelList(:,3))*INFO.SliceThickness,30,'ro')
        for i = 1:length(Matched)
            if length(Matched{i})>0
                closeMatch{i} = sqrt((mean(data_subpl(j).PixelList(:,1))-Matched{i}(:,1)).^2+(mean(data_subpl(j).PixelList(:,2))-Matched{i}(:,2)).^2+(mean(data_subpl(j).PixelList(:,3))-Matched{i}(:,3)).^2);
                lengthMatch{i} = [0;cumsum(abs(closeMatch{i}(1:(end-1)) - closeMatch{i}(2:end)))];
            else
                closeMatch{i} = sqrt(sum((size(Images).*[INFO.PixelSpacing(1) INFO.PixelSpacing(2) INFO.SliceThickness]).^2));
                lengthMatch{i} = [0;cumsum(abs(closeMatch{i}(1:(end-1)) - closeMatch{i}(2:end)))];
            end
        end
        whichClose = (cell2mat(cellfun(@min,closeMatch,'uni',false))<(RRadMult)).*(1:length(Matched));
        whichClose(whichClose==0) = [];
        for i = 1:length(whichClose)
            [~,Pos] = min(closeMatch{whichClose(i)});
            notClose = 1:length(Matched);
            notClose(whichClose(i)) = [];
            OMatched = vertcat(Matched{notClose});
            CMatched = Matched{whichClose(i)};
            scatter3(CMatched(Pos,1)*INFO.PixelSpacing(1),CMatched(Pos,2)*INFO.PixelSpacing(2),CMatched(Pos,3)*INFO.SliceThickness,15,'ko','filled')
            DeleteClearance = 0;
            if abs(lengthMatch{whichClose(i)}(end)-lengthMatch{whichClose(i)}(Pos))<Dclose
                Pos = length(closeMatch{whichClose(i)});
                closeMatch{whichClose(i)} = flipud([0;cumsum(abs(flipud(closeMatch{whichClose(i)}(2:(end)))-flipud(closeMatch{whichClose(i)}(1:(end-1)))))]);
                if (sum(sqrt(sum(((CMatched(Pos,1:3)-OMatched(:,1:3)).*[INFO.PixelSpacing' INFO.SliceThickness]).^2,2))<SmallLinks)==0)>0
                    DeleteClearance = 1;
                end
            elseif abs(lengthMatch{whichClose(i)}(1)-lengthMatch{whichClose(i)}(Pos))<Dclose
                Pos = 1;
                closeMatch{whichClose(i)} = [0;cumsum(abs(closeMatch{whichClose(i)}(1:(end-1))-closeMatch{whichClose(i)}(2:(end))))];
                if (sum(sqrt(sum(((CMatched(Pos,1:3)-OMatched(:,1:3)).*[INFO.PixelSpacing' INFO.SliceThickness]).^2,2))<SmallLinks)==0)>0
                    abs(closeMatch{whichClose(i)}(1)-closeMatch{whichClose(i)}(Pos))
                    DeleteClearance = 1;
                end
            end
            if DeleteClearance>0
                [~,Idx] = sort(closeMatch{whichClose(i)});
                try
                    for k = 1:length(Idx)
                        if (sum(sqrt(sum(((CMatched(Idx(k),1:3)-OMatched(:,1:3)).*[INFO.PixelSpacing' INFO.SliceThickness]).^2,2))<SmallLinks)+...
                                (CMatched(Idx(k),3)<=BottomBranch))>0
                            break
                        end
                        scatter3(Matched{whichClose(i)}(Idx(k),1)*INFO.PixelSpacing(1),Matched{whichClose(i)}(Idx(k),2)*INFO.PixelSpacing(2),Matched{whichClose(i)}(Idx(k),3)*INFO.SliceThickness,10,'ro','filled')
                        Matched{whichClose(i)}(Idx(k),:) = [];
                    end
                end
            end
        end
    end
    pause(0.01)
    FullList = vertcat(Matched{:});
    blanki3 = zeros(size(blanki));
    blanki3(sub2ind(size(Images),FullList(:,2),FullList(:,1),FullList(:,3))) = 1;
    blanki3 = imdilate(blanki3,strel('sphere',round(1)));
    blanksr = bwconncomp(blanki3,6);
    data_full = regionprops(blanksr, 'PixelList');
    VoxelList = cell2mat(cellfun(@length,cellfun(@transpose,{data_full.PixelList},'UniformOutput',false),'uni',false));
    if sum(VoxelList<(max(VoxelList)*0.01))>0
        RedivideList = ((VoxelList<(max(VoxelList)*0.01)).*(1:length(VoxelList)));
        RedivideList(RedivideList==0) = [];
        for j = RedivideList
            FullList(ismember(sub2ind(size(blanki),FullList(:,1),FullList(:,2),FullList(:,3)),...
                sub2ind(size(blanki),data_full(j).PixelList(:,1),data_full(j).PixelList(:,2),...
                data_full(j).PixelList(:,3))),:) = [];
        end
    end
    %%the end of opional mask refinement/culling
    
    %%the outlier detection needs to be redone if culling was performed. 
    blanksr = bwconncomp(blanki,6);
    data_full = regionprops(blanksr, 'PixelList');
    data_full = vertcat(data_full.PixelList);
    data_full2 = zeros(length(data_full),4);
    for i=1:length(data_full)
        dist2 = sum((FullList(:,1:3) - data_full(i,:)).^2,2);
        data_full2(i,:) = [data_full(i,:),min(FullList(dist2 == min(dist2),4))];
    end
    blanks = zeros(size(blanki));
    blanks(sub2ind(size(Images),FullList(:,2),FullList(:,1),FullList(:,3))) = 1;
    blankp = zeros(size(blanki));
    blankp(sub2ind(size(Images),data_full2(:,2),data_full2(:,1),data_full2(:,3))) = data_full2(:,4);
    blankd = bwdistsc(imcomplement(blanki),[INFO.PixelSpacing(1) INFO.PixelSpacing(2) INFO.SliceThickness]);
    blankdiam = bwdistsc(blanks,[INFO.PixelSpacing(1) INFO.PixelSpacing(2) INFO.SliceThickness]).*blanki;
    
    data_pl = data_pl(logical((data_pl(:,3)>size(Images,3)*0.1).*(data_pl(:,3)<size(Images,3)*0.95)),:);
    [fitresult1] = fit([blankp(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3))),...
        blankd(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3)))],...
        blankdiam(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3))),'poly41');
    
    [pci py] = predint(fitresult1,[blankp(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3))) ...
        blankd(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3)))],ConInt2,'observation','off');
    ry = blankdiam(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3)));
    SubSec = data_pl(logical(ry>pci(:,2)),:);
    
    blankSub = zeros(size(Images));
    blankSub(sub2ind(size(Images),SubSec(:,2), SubSec(:,1), SubSec(:,3))) = 1;
    data_subc = bwconncomp(blankSub,6);
    data_subpl = regionprops(data_subc, 'PixelList');
    Voxels_sub = cell2mat(cellfun(@length,cellfun(@transpose,{data_subpl.PixelList},'UniformOutput',false),'uni',false));
    
    try
        i = 1;
        while i <= size(data_subpl,1)
            OutlierProps = mean(fitresult1(blankp(sub2ind(size(Images),data_subpl(i).PixelList(:,2),data_subpl(i).PixelList(:,1),data_subpl(i).PixelList(:,3))),...
                blankd(sub2ind(size(Images),data_subpl(i).PixelList(:,2),data_subpl(i).PixelList(:,1),data_subpl(i).PixelList(:,3)))));
            OutlierProps(OutlierProps<0) = 0;
            if Voxels_sub(i)<(OutlierProps*minClustRate(1)+minClustRate(2))
                data_subpl(i) = [];
                Voxels_sub(i) = [];
            else
                i = i+1;
            end
        end
    end
    %%end of repeat outlier analysis when culling occurs. 
end

outlierProps = cell(1,1);
for i = 1:length(data_subpl)
    dList = blankd(sub2ind(size(Images),data_subpl(i).PixelList(:,2),data_subpl(i).PixelList(:,1),data_subpl(i).PixelList(:,3)));
    diamList = blankdiam(sub2ind(size(Images),data_subpl(i).PixelList(:,2),data_subpl(i).PixelList(:,1),data_subpl(i).PixelList(:,3)));
    hlist = blankh(sub2ind(size(Images),data_subpl(i).PixelList(:,2),data_subpl(i).PixelList(:,1),data_subpl(i).PixelList(:,3)));
    pList = blankp(sub2ind(size(Images),data_subpl(i).PixelList(:,2),data_subpl(i).PixelList(:,1),data_subpl(i).PixelList(:,3)));
    eList = blankdiam(sub2ind(size(Images),data_subpl(i).PixelList(:,2),data_subpl(i).PixelList(:,1),data_subpl(i).PixelList(:,3)))-...
        fitresult1(blankp(sub2ind(size(Images),data_subpl(i).PixelList(:,2),data_subpl(i).PixelList(:,1),data_subpl(i).PixelList(:,3))),...
        blankd(sub2ind(size(Images),data_subpl(i).PixelList(:,2),data_subpl(i).PixelList(:,1),data_subpl(i).PixelList(:,3))));
    outlierProps{i} = [hlist,dList,diamList,pList,eList];
end
SubList = vertcat(data_subpl.PixelList);


figure
minH = min(blankh(blankh>0));
scatter3(blankp(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3))),...
    blankd(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3))),...
    blankdiam(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3))),10,...
    blankh(sub2ind(size(Images),data_pl(:,2),data_pl(:,1),data_pl(:,3))) - minH,'filled',...
    'MarkerFaceAlpha', 0.3,'MarkerEdgeAlpha', 0.3)
view(135,30);
xlabel('Base Distance {\it d_b} (mm)')
ylabel('Edge Distance {\it d_e} (mm)')
zlabel('Centreline Distance {\it d_c} (mm)')
set(gca, 'YAxisLocation', 'right')
c = colorbar('location','EastOutside');
c.Label.String = 'Time of Flight (value)';
hold on
xlims = get(gca,'XLim');
ylims = get(gca,'YLim');
zlims = get(gca,'ZLim');
x = repmat(linspace(xlims(1),xlims(2),50),50,1);
y = repmat(linspace(ylims(1),ylims(2),50)',1,50);
surf(x,y,fitresult1(x,y),'FaceAlpha',0)
xlim([0 xlims(2)]) %325
ylim([0 ylims(2)])
zlim([0 zlims(2)]) %8
c.TickLabels = cellstr(num2str(str2double(c.TickLabels)+minH));


if ~isempty(SubList)
    scatter3(blankp(sub2ind(size(Images),SubList(:,2),SubList(:,1),SubList(:,3))),...
        blankd(sub2ind(size(Images),SubList(:,2),SubList(:,1),SubList(:,3))),...
        blankdiam(sub2ind(size(Images),SubList(:,2),SubList(:,1),SubList(:,3))),10,'ro','filled')
end

UIAsize = Voxels_sub*prod([INFO.PixelSpacing' INFO.SliceThickness]);
for i=1:size(UIAsize,2)
    strcat('UIA discovered! Size: ',num2str(UIAsize(i)),{' '},'mm3', '. Dim:',{' '},num2str([max(data_subpl(i).PixelList)-min(data_subpl(i).PixelList)].*[INFO.PixelSpacing' INFO.SliceThickness]),{' '},'mm.')
end
end

