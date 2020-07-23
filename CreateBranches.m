function [SubInd,blanks,blanki] = CreateBranches(Images,blanks,blankd,blanki,INFO,SmallLinks,NumTissue,WhichTissue,Kw,binnum,ObjectImageFrac,MMThresh,RegionRad,PeriphPrc,TightenSpacing,reThresh)

dist = 0;
blankir = bwconncomp(blanki,6);
data_plall1 = regionprops(blankir, 'PixelList');
data_plall1 = vertcat(data_plall1.PixelList);
blanksr = bwconncomp(blanks,6);
data_pls = regionprops(blanksr, 'PixelList');
data_pls = vertcat(data_pls.PixelList);
minZ = 1;
try
    while ~ismember(1,blanks(:,:,minZ))
        minZ = minZ + 1;
    end
end
[~,Loc] = ismember(1,blanks(:,:,minZ));
[minX,minY] = ind2sub(size(blanks(:,:,minZ)),Loc);
SubInd{1} = [minY,minX,minZ,dist];

figure
scatter3(data_plall1(:,1)*INFO.PixelSpacing(1),data_plall1(:,2)*INFO.PixelSpacing(2),data_plall1(:,3)*INFO.SliceThickness,10,'k','filled',...
    'MarkerEdgeAlpha', 0.02,'MarkerFaceAlpha', 0.02)
view(180,90);
hold on
scatter3(data_pls(:,1)*INFO.PixelSpacing(1),data_pls(:,2)*INFO.PixelSpacing(2),data_pls(:,3)*INFO.SliceThickness,'k.')
    xlabel('mm')
    ylabel('mm')
    zlabel('mm')
pause(0.01)

ThreshMap = zeros(size(Images));

k=1;
j=1;
try
    while size(data_pls,1)>0
        lastThresh = 1;
        while dist<SmallLinks
            data_pls(ismember(sub2ind(size(Images),data_pls(:,1),data_pls(:,2),data_pls(:,3)),...
                sub2ind(size(Images),SubInd{k}(j,1),SubInd{k}(j,2),SubInd{k}(j,3))),:) = [];
            [dist,which] = min(sqrt(((SubInd{k}(j,1)-data_pls(:,1)).*INFO.PixelSpacing(1)).^2+...
                ((SubInd{k}(j,2)-data_pls(:,2)).*INFO.PixelSpacing(2)).^2+...
                ((SubInd{k}(j,3)-data_pls(:,3)).*INFO.SliceThickness).^2));
            NumAdded = 0;
            j = j+NumAdded+1;
            SubInd{k}(j,:) = [data_pls(which,:),SubInd{k}(end-NumAdded,4)+dist];
            MaxSpacing = sqrt(4*RegionRad^2-(blankd(SubInd{k}(j,2),SubInd{k}(j,1),SubInd{k}(j,3)))^2); %calculating d from 2a and R in: http://mathworld.wolfram.com/Sphere-SphereIntersection.html
            if ((sqrt(((SubInd{k}(j-1,1)-SubInd{k}(lastThresh,1)).*INFO.PixelSpacing(1)).^2+...
                    ((SubInd{k}(j-1,2)-SubInd{k}(lastThresh,2)).*INFO.PixelSpacing(2)).^2+...
                    ((SubInd{k}(j-1,3)-SubInd{k}(lastThresh,3)).*INFO.SliceThickness).^2)>(MaxSpacing*TightenSpacing))+(lastThresh==1))*(reThresh==1)>0
                
                LocalRegion = zeros(size(Images));
                LocalRegion(sub2ind(size(Images),SubInd{k}(end-1,2),SubInd{k}(end-1,1),SubInd{k}(end-1,3))) = 1;
                LocalRegion = bwdistsc(LocalRegion,[INFO.PixelSpacing', INFO.SliceThickness])<RegionRad;
                RegionVox = find(LocalRegion>0);
                LocalRegion = int16(LocalRegion).*int16(Images);
                LocalRegion = double(LocalRegion);
                LocalRegion(LocalRegion==0) = NaN;
               
                [Threshold,~,~] = AUTOthresh(NumTissue,WhichTissue,Kw,LocalRegion,binnum,ObjectImageFrac,PeriphPrc,MMThresh,0);
                
                blanklr = bwconncomp(~isnan(LocalRegion),6);
                data_lr = regionprops(blanklr, 'PixelList');
                data_lr = vertcat(data_lr.PixelList);
                 scatter3(SubInd{k}(j-1,1)*INFO.PixelSpacing(1),SubInd{k}(j-1,2)*INFO.PixelSpacing(2),SubInd{k}(j-1,3)*INFO.SliceThickness,RegionRad*(150*INFO.PixelSpacing(1)/0.38),Threshold(1),...
                    'filled','MarkerEdgeAlpha', 0.1,'MarkerFaceAlpha', 0.1)
                text(SubInd{k}(j-1,1)*INFO.PixelSpacing(1),SubInd{k}(j-1,2)*INFO.PixelSpacing(2),SubInd{k}(j-1,3)*INFO.SliceThickness,num2str(Threshold(1),'%4.0f'),...
                    'color',[0.5 0.5 0.5],'FontSize',6,'HorizontalAlignment', 'Center')
                pause(0.01)
                
                [LocalRegion,~] = DICOMthresh(LocalRegion,Threshold(1),0,0);
                ThreshMap(RegionVox) = Threshold(1);
                blanki(RegionVox) = LocalRegion(RegionVox);
                lastThresh = j-1;
            end
        end
        
        
        if (((j-1)~=lastThresh)*(reThresh==1))>0
            LocalRegion = zeros(size(Images));
            LocalRegion(sub2ind(size(Images),SubInd{k}(end-1,2),SubInd{k}(end-1,1),SubInd{k}(end-1,3))) = 1;
            LocalRegion = bwdistsc(LocalRegion,[INFO.PixelSpacing', INFO.SliceThickness])<RegionRad;
            RegionVox = find(LocalRegion>0);
            LocalRegion = int16(LocalRegion).*int16(Images);
            LocalRegion = double(LocalRegion);
            LocalRegion(LocalRegion==0) = NaN;
            [Threshold,~,~] = AUTOthresh(NumTissue,WhichTissue,Kw,LocalRegion,binnum,ObjectImageFrac,PeriphPrc,MMThresh,0);
            blanklr = bwconncomp(~isnan(LocalRegion),6);
            data_lr = regionprops(blanklr, 'PixelList');
            data_lr = vertcat(data_lr.PixelList);
            scatter3(SubInd{k}(j-1,1)*INFO.PixelSpacing(1),SubInd{k}(j-1,2)*INFO.PixelSpacing(2),SubInd{k}(j-1,3)*INFO.SliceThickness,RegionRad*(150*INFO.PixelSpacing(1)/0.38),Threshold(1),...
                'filled','MarkerEdgeAlpha', 0.1,'MarkerFaceAlpha', 0.1)
            text(SubInd{k}(j-1,1)*INFO.PixelSpacing(1),SubInd{k}(j-1,2)*INFO.PixelSpacing(2),SubInd{k}(j-1,3)*INFO.SliceThickness,num2str(Threshold(1),'%4.0f'),...
                'color',[0.5 0.5 0.5],'FontSize',6,'HorizontalAlignment', 'Center')
            pause(0.01)
            [LocalRegion,~] = DICOMthresh(LocalRegion,Threshold(1),0,0);
            ThreshMap(RegionVox) = Threshold(1);
            blanki(RegionVox) = LocalRegion(RegionVox);
        end
        SubInd{k}(end,:) = [];
        
        k = k+1;
        totalRecord = vertcat(SubInd{:});
        totalDIST = sqrt(((repmat(totalRecord(:,1)',size(data_pls,1),1) - repmat(data_pls(:,1),1,size(totalRecord,1))).*INFO.PixelSpacing(1)).^2+...
            ((repmat(totalRecord(:,2)',size(data_pls,1),1) - repmat(data_pls(:,2),1,size(totalRecord,1))).*INFO.PixelSpacing(2)).^2+...
            ((repmat(totalRecord(:,3)',size(data_pls,1),1) - repmat(data_pls(:,3),1,size(totalRecord,1))).*INFO.SliceThickness).^2);
        [dist,linearIND] = min(totalDIST(:));
        [plsNum,SubNum] = ind2sub(size(totalDIST,1),linearIND);
        for l = 1:length(SubInd)
            SubNum = SubNum - size(SubInd{l},1);
            if SubNum<=0
                SubNum = SubNum + size(SubInd{l},1);
                SubCell = l;
                break
            end
        end
        SubInd{k}(1:2,:) = [SubInd{SubCell}(SubNum,1:3),0;...
            data_pls(plsNum,:),dist];
        j = 2;
        if dist>SmallLinks 
            data_pls(ismember(sub2ind(size(Images),data_pls(:,1),data_pls(:,2),data_pls(:,3)),...
                sub2ind(size(Images),SubInd{k}(1,1),SubInd{k}(1,2),SubInd{k}(1,3))),:) = []; %this is the deleter
            dist = 0;
            SubInd{k}(1,:) = [];
            SubInd{k}(1,4) = dist;
            j = 1;
        end
    end
end
xlabel('mm')
ylabel('mm')
zlabel('mm')
view(180,90)
colormap(jet)


blankir = bwconncomp(blanki,6);
data_plall = regionprops(blankir, 'PixelList');
Voxels_sub = cell2mat(cellfun(@length,cellfun(@transpose,{data_plall.PixelList},'UniformOutput',false),'uni',false));
Voxels_sub = (Voxels_sub>(numel(Images)*ObjectImageFrac)).*(1:length(Voxels_sub)); %ObjectImageFrac used to be reduced by 0.1, needed 0.25 for ABN113, but 1 (normal) should work.
Voxels_sub(Voxels_sub==0) = [];
data_plall = vertcat(data_plall(Voxels_sub).PixelList);

%% Global Rethresh
if reThresh>0
    blanke = (ThreshMap==0);
    [~,idx] = bwdist(ThreshMap);
    blanke = ThreshMap(idx).*blanke;
    blanke(blanke==0) = NaN;
    Thresh1 = (Images>blanke);
    Thresh1(isnan(Thresh1)) = 0;
    blanki2 = Thresh1 + blanki;
    blankir2 = bwconncomp(blanki2,6);
    data_plall2 = regionprops(blankir2, 'PixelList');
    for i = 1:length(data_plall2)
        if sum(ismember(sub2ind(size(Images),data_plall(:,2),data_plall(:,1),data_plall(:,3)),...
                sub2ind(size(Images),data_plall2(i).PixelList(:,2),data_plall2(i).PixelList(:,1),data_plall2(i).PixelList(:,3))))>0
            data_plall(ismember(sub2ind(size(Images),data_plall(:,2),data_plall(:,1),data_plall(:,3)),...
                sub2ind(size(Images),data_plall2(i).PixelList(:,2),data_plall2(i).PixelList(:,1),data_plall2(i).PixelList(:,3))),:) = [];
            data_plall = vertcat(data_plall,data_plall2(i).PixelList);
        end
    end
    
    %% redefine and resegment all this ish
    blanki2 = zeros(size(Images));
    blanki2(sub2ind(size(Images),data_plall(:,2),data_plall(:,1),data_plall(:,3))) = 1;
    scatter3(data_plall(:,1)*INFO.PixelSpacing(1),data_plall(:,2)*INFO.PixelSpacing(2),data_plall(:,3)*INFO.SliceThickness,10,'k','filled',...
        'MarkerEdgeAlpha', 0.1,'MarkerFaceAlpha', 0.1)
    blanks = bwskel(imclose(imdilate(bwskel(imfill(blanki2>0,4,'holes')),strel('sphere',1)),strel('sphere',1))); %dilating at 2 looks cleaner - but has downstream problems
    blanksr = bwconncomp(blanks,6);
    data_sl = regionprops(blanksr, 'PixelList');
    data_sl = vertcat(data_sl.PixelList);
    scatter3(data_sl(:,1)*INFO.PixelSpacing(1),data_sl(:,2)*INFO.PixelSpacing(2),data_sl(:,3)*INFO.SliceThickness,'r.')
    scatter3(data_pls(:,1)*INFO.PixelSpacing(1),data_pls(:,2)*INFO.PixelSpacing(2),data_pls(:,3)*INFO.SliceThickness,'k.')
    scatter3(data_plall(:,1)*INFO.PixelSpacing(1),data_plall(:,2)*INFO.PixelSpacing(2),data_plall(:,3)*INFO.SliceThickness,10,'k','filled',...
    'MarkerEdgeAlpha', 0.05,'MarkerFaceAlpha', 0.05)
    
    blanki = blanki2;
end

end