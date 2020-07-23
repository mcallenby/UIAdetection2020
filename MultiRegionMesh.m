function MultiRegionMesh(subRegions,Images,Threshold,INFO,data_pl,blankp,blankd,blankdiam,blankh,data_plsi,data_skel,colorstring,ChoiceDisplay, fitresult)

UniqueSubRegions = [];
for i = 1:size(subRegions,1)
    for j = 1:size(subRegions,2)
        if length(subRegions{i,j})>0
            UniqueSubRegions(size(UniqueSubRegions,1)+1,:) = subRegions{i,j};
            for l = 1:size(UniqueSubRegions,1)
                ExpVec(l) = ~isempty(intersect(UniqueSubRegions(end,1):UniqueSubRegions(end,2),UniqueSubRegions(l,1):UniqueSubRegions(l,2)))*...
                    ~isempty(intersect(UniqueSubRegions(end,3):UniqueSubRegions(end,4),UniqueSubRegions(l,3):UniqueSubRegions(l,4)))*...
                    ~isempty(intersect(UniqueSubRegions(end,5):UniqueSubRegions(end,6),UniqueSubRegions(l,5):UniqueSubRegions(l,6)));
            end
            ExpVec = ExpVec.*(1:length(ExpVec));
            ExpVec(ExpVec == 0) = [];
            ExpVec(ExpVec == size(UniqueSubRegions,1)) = [];
            if ~isempty(ExpVec)
                for l = 1:6
                    Expansion = UniqueSubRegions(end,l) - UniqueSubRegions(ExpVec,l);
                    if (min(Expansion)<0) && (l == 1 || l == 3 || l == 5)
                        UniqueSubRegions(ExpVec,l) = UniqueSubRegions(ExpVec,l) + min(Expansion);
                    elseif (max(Expansion)>0) && (l == 2 || l == 4 || l == 6)
                        UniqueSubRegions(ExpVec,l) = UniqueSubRegions(ExpVec,l) + max(Expansion);
                    end
                end
                UniqueSubRegions(end,:) = [];
            end
        end
    end
end

ImageIndex = 1:length(Images);
for i = ImageIndex
    ImageHeight(i) = size(Images{i},3);
end
[~,tallImage] = max(ImageHeight);
ImageIndex(ImageIndex==tallImage) = [];

%plot a timecourse of the troublesome areas
[optimizer, metric] = imregconfig('monomodal');
for i=1:size(UniqueSubRegions,1)
    figure
    hold on
    for j=[tallImage,ImageIndex]
        ImagesSub2{j} = zeros(size(Images{j}));
        
        xMat{i} = round(UniqueSubRegions(i,3):UniqueSubRegions(i,4));
        yMat{i} = round(UniqueSubRegions(i,1):UniqueSubRegions(i,2));
        zMat{i} = round(UniqueSubRegions(i,5):UniqueSubRegions(i,6));
        xMat{i}(((xMat{i} <= 0) + (xMat{i} >= size(Images{j},2)))>0) = [];
        yMat{i}(((yMat{i} <= 0) + (yMat{i} >= size(Images{j},1)))>0) = [];
        zMat{i}(((zMat{i} <= 0) + (zMat{i} >= size(Images{j},3)))>0) = [];
        
        ImagesSub{j} = Images{j}(xMat{i},yMat{i},zMat{i});
        if ismember(j,ImageIndex)
            tform = imregtform(ImagesSub{j},ImagesSub{tallImage},'affine',optimizer, metric);
            ImagesSub{j} = imwarp(ImagesSub{j},tform,'OutputView',imref3d(size(ImagesSub{j})));
        end
        ImagesSub2{j}(xMat{i},yMat{i},zMat{i}) = ImagesSub{j};
        [~,data_plj] = DICOMthresh(ImagesSub2{j},Threshold{j},0,0);
        if (~isempty(data_plj))*(j==1)>0
            scatter3(data_plj(:,1)*INFO.PixelSpacing(1),data_plj(:,2)*INFO.PixelSpacing(2),data_plj(:,3)*INFO.SliceThickness,180,...
                blankh{j}(sub2ind(size(Images{j}),data_plj(:,2),data_plj(:,1),data_plj(:,3))),'.','MarkerEdgeAlpha', 0.1,'MarkerFaceAlpha', 0.1)
        elseif ~isempty(data_plj)
            scatter3(data_plj(:,1)*INFO.PixelSpacing(1),data_plj(:,2)*INFO.PixelSpacing(2),data_plj(:,3)*INFO.SliceThickness,180,colorstring{j},...
                'MarkerEdgeAlpha', 0.1,'MarkerFaceAlpha', 0.1)
        end
        scatter3(data_plsi{ChoiceDisplay}(:,1)*INFO.PixelSpacing(2),data_plsi{ChoiceDisplay}(:,2)*INFO.PixelSpacing(1),...
            data_plsi{ChoiceDisplay}(:,3)*INFO.SliceThickness,180,'r.','MarkerEdgeAlpha', 0.5,'MarkerFaceAlpha', 0.5)
        axis([UniqueSubRegions(i,1)*INFO.PixelSpacing(2) UniqueSubRegions(i,2)*INFO.PixelSpacing(2) UniqueSubRegions(i,3)*INFO.PixelSpacing(1) UniqueSubRegions(i,4)*INFO.PixelSpacing(1) UniqueSubRegions(i,5)*INFO.SliceThickness UniqueSubRegions(i,6)*INFO.SliceThickness])
    end
    xlabel('mm')
    ylabel('mm')
    zlabel('mm')
end


%plot the final scan and highlight troublesome areas to track over time
for i = 1:length(ChoiceDisplay)
    
    figure
    scatter3(data_pl{ChoiceDisplay(i)}(:,1)*INFO.PixelSpacing(1),data_pl{ChoiceDisplay(i)}(:,2)*INFO.PixelSpacing(2),data_pl{ChoiceDisplay(i)}(:,3)*INFO.SliceThickness,10,'k','filled')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.45]);
    view(180,90);
    title(strcat(INFO.ContentDate,{': '},INFO.StudyDescription))
    xlabel('mm')
    ylabel('mm')
    zlabel('mm')
    set(gca, 'YAxisLocation', 'right')
    
            figure
    scatter3(data_pl{ChoiceDisplay(i)}(:,1)*INFO.PixelSpacing(1),data_pl{ChoiceDisplay(i)}(:,2)*INFO.PixelSpacing(2),data_pl{ChoiceDisplay(i)}(:,3)*INFO.SliceThickness,10,...
        blankh{ChoiceDisplay(i)}(sub2ind(size(Images{ChoiceDisplay(i)}),data_pl{ChoiceDisplay(i)}(:,2),data_pl{ChoiceDisplay(i)}(:,1),data_pl{ChoiceDisplay(i)}(:,3))),'filled',...
        'MarkerEdgeAlpha', 0.2,'MarkerFaceAlpha', 0.2)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.45]);
    view(180,90);
    title(strcat(INFO.ContentDate,{': '},INFO.StudyDescription))
    xlabel('mm')
    ylabel('mm')
    zlabel('mm')
    set(gca, 'YAxisLocation', 'right')
    c = colorbar('location','WestOutside');
    c.Label.String = 'Time of Flight (units)';
    
    figure
    scatter3(data_pl{ChoiceDisplay(i)}(:,1)*INFO.PixelSpacing(1),data_pl{ChoiceDisplay(i)}(:,2)*INFO.PixelSpacing(2),data_pl{ChoiceDisplay(i)}(:,3)*INFO.SliceThickness,10,...
        blankd{ChoiceDisplay(i)}(sub2ind(size(Images{ChoiceDisplay(i)}),data_pl{ChoiceDisplay(i)}(:,2),data_pl{ChoiceDisplay(i)}(:,1),data_pl{ChoiceDisplay(i)}(:,3))),'filled',...
        'MarkerEdgeAlpha', 0.2,'MarkerFaceAlpha', 0.2)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.45]);
    view(180,90);
    title(strcat(INFO.ContentDate,{': '},INFO.StudyDescription))
    xlabel('mm')
    ylabel('mm')
    zlabel('mm')
    set(gca, 'YAxisLocation', 'right')
    c = colorbar('location','WestOutside');
    c.Label.String = 'Edge Distance {\it d_e} (mm)';
    
        figure
    scatter3(data_pl{ChoiceDisplay(i)}(:,1)*INFO.PixelSpacing(1),data_pl{ChoiceDisplay(i)}(:,2)*INFO.PixelSpacing(2),data_pl{ChoiceDisplay(i)}(:,3)*INFO.SliceThickness,10,...
        blankdiam{ChoiceDisplay(i)}(sub2ind(size(Images{ChoiceDisplay(i)}),data_pl{ChoiceDisplay(i)}(:,2),data_pl{ChoiceDisplay(i)}(:,1),data_pl{ChoiceDisplay(i)}(:,3))),'filled',...
        'MarkerEdgeAlpha', 0.1,'MarkerFaceAlpha', 0.1)
    %hold on
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.45]);
    view(180,90);
    title(strcat(INFO.ContentDate,{': '},INFO.StudyDescription))
    xlabel('mm')
    ylabel('mm')
    zlabel('mm')
    set(gca, 'YAxisLocation', 'right')
    c = colorbar('location','WestOutside');
    c.Label.String = 'Centreline Distance {\it d_c} (mm)';
    %scatter3(data_skel{ChoiceDisplay(i)}(:,1)*INFO.PixelSpacing(1),data_skel{ChoiceDisplay(i)}(:,2)*INFO.PixelSpacing(2),data_skel{ChoiceDisplay(i)}(:,3)*INFO.SliceThickness,5,'k','filled')
    
            figure
    scatter3(data_pl{ChoiceDisplay(i)}(:,1)*INFO.PixelSpacing(1),data_pl{ChoiceDisplay(i)}(:,2)*INFO.PixelSpacing(2),data_pl{ChoiceDisplay(i)}(:,3)*INFO.SliceThickness,10,...
        blankp{ChoiceDisplay(i)}(sub2ind(size(Images{ChoiceDisplay(i)}),data_pl{ChoiceDisplay(i)}(:,2),data_pl{ChoiceDisplay(i)}(:,1),data_pl{ChoiceDisplay(i)}(:,3))),'filled',...
        'MarkerEdgeAlpha', 0.5,'MarkerFaceAlpha', 0.5)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.45]);
    view(180,90);
    title(strcat(INFO.ContentDate,{': '},INFO.StudyDescription))
    xlabel('mm')
    ylabel('mm')
    zlabel('mm')
    set(gca, 'YAxisLocation', 'right')
    c = colorbar('location','WestOutside');
    c.Label.String = 'Base Distance {\it d_b} (mm)';
    
    figure
    scatter3(data_pl{ChoiceDisplay(i)}(:,1)*INFO.PixelSpacing(1),data_pl{ChoiceDisplay(i)}(:,2)*INFO.PixelSpacing(2),data_pl{ChoiceDisplay(i)}(:,3)*INFO.SliceThickness,10,...
    (blankdiam{ChoiceDisplay(i)}(sub2ind(size(Images{ChoiceDisplay(i)}),data_pl{ChoiceDisplay(i)}(:,2),data_pl{ChoiceDisplay(i)}(:,1),data_pl{ChoiceDisplay(i)}(:,3)))-...
    fitresult(blankp{ChoiceDisplay(i)}(sub2ind(size(Images{ChoiceDisplay(i)}),data_pl{ChoiceDisplay(i)}(:,2),data_pl{ChoiceDisplay(i)}(:,1),data_pl{ChoiceDisplay(i)}(:,3))),...
    blankd{ChoiceDisplay(i)}(sub2ind(size(Images{ChoiceDisplay(i)}),data_pl{ChoiceDisplay(i)}(:,2),data_pl{ChoiceDisplay(i)}(:,1),data_pl{ChoiceDisplay(i)}(:,3))))),...
    'filled','MarkerEdgeAlpha', 0.5,'MarkerFaceAlpha', 0.5)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.45]);
    view(180,90);
    title(strcat(INFO.ContentDate,{': '},INFO.StudyDescription))
    xlabel('mm')
    ylabel('mm')
    zlabel('mm')
    set(gca, 'YAxisLocation', 'right')
    c = colorbar('location','WestOutside');
    caxis([ 0 , 3])
    c.Label.String = 'Regression Error (mm)';
    
    
        figure
    scatter3(data_pl{ChoiceDisplay(i)}(:,1)*INFO.PixelSpacing(1),data_pl{ChoiceDisplay(i)}(:,2)*INFO.PixelSpacing(2),data_pl{ChoiceDisplay(i)}(:,3)*INFO.SliceThickness,10,...
        blankh{ChoiceDisplay(i)}(sub2ind(size(Images{ChoiceDisplay(i)}),data_pl{ChoiceDisplay(i)}(:,2),data_pl{ChoiceDisplay(i)}(:,1),data_pl{ChoiceDisplay(i)}(:,3))),'filled',...
        'MarkerEdgeAlpha', 0.1,'MarkerFaceAlpha', 0.1)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.45]*1.45);
    view(180,90);
    title(strcat(INFO.ContentDate,{': '},INFO.StudyDescription))
    xlabel('mm')
    ylabel('mm')
    zlabel('mm')
    set(gca, 'YAxisLocation', 'right')
    c = colorbar('location','WestOutside');
    c.Label.String = 'Time of Flight (units)';
    hold on
    
    %scatter3 the skeleton on this as well
    if ~isempty(data_plsi{ChoiceDisplay(i)})
       %%important if you are lining up many choices, otherwise no worries
%         index = convhull(data_plsi{ChoiceDisplay(i)}(:,1)*INFO.PixelSpacing(1),data_plsi{ChoiceDisplay(i)}(:,2)*INFO.PixelSpacing(2),...
%             data_plsi{ChoiceDisplay(i)}(:,3)*INFO.SliceThickness);
%         scatter3(data_plsi{ChoiceDisplay(i)}(index,1)*INFO.PixelSpacing(1),data_plsi{ChoiceDisplay(i)}(index,2)*INFO.PixelSpacing(2),...
%             data_plsi{ChoiceDisplay(i)}(index,3)*INFO.SliceThickness,'r.','MarkerEdgeAlpha', 0.5,'MarkerFaceAlpha', 0.5)
                scatter3(data_plsi{ChoiceDisplay(i)}(:,1)*INFO.PixelSpacing(1),data_plsi{ChoiceDisplay(i)}(:,2)*INFO.PixelSpacing(2),...
            data_plsi{ChoiceDisplay(i)}(:,3)*INFO.SliceThickness,'r.','MarkerEdgeAlpha', 0.5,'MarkerFaceAlpha', 0.5)
    end
    
    for i=1:size(UniqueSubRegions,1)
        X1 = [UniqueSubRegions(i,1);UniqueSubRegions(i,1);UniqueSubRegions(i,2);UniqueSubRegions(i,2);...
            UniqueSubRegions(i,2);UniqueSubRegions(i,2);UniqueSubRegions(i,1);UniqueSubRegions(i,1);...
            UniqueSubRegions(i,1);UniqueSubRegions(i,1)];
        Y1 = [UniqueSubRegions(i,3);UniqueSubRegions(i,3);UniqueSubRegions(i,3);UniqueSubRegions(i,3);...
            UniqueSubRegions(i,4);UniqueSubRegions(i,4);UniqueSubRegions(i,4);UniqueSubRegions(i,4);...
            UniqueSubRegions(i,3);UniqueSubRegions(i,3)];
        Z1 = [UniqueSubRegions(i,5);UniqueSubRegions(i,6);UniqueSubRegions(i,6);UniqueSubRegions(i,5);...
            UniqueSubRegions(i,5);UniqueSubRegions(i,6);UniqueSubRegions(i,6);UniqueSubRegions(i,5);...
            UniqueSubRegions(i,5);UniqueSubRegions(i,6)];
        Z2 = [UniqueSubRegions(i,6);UniqueSubRegions(i,5);UniqueSubRegions(i,5);UniqueSubRegions(i,6);...
            UniqueSubRegions(i,6);UniqueSubRegions(i,5);UniqueSubRegions(i,5);UniqueSubRegions(i,6);...
            UniqueSubRegions(i,6);UniqueSubRegions(i,5)];
        
        plot3(X1*INFO.PixelSpacing(1),Y1*INFO.PixelSpacing(2),Z1*INFO.SliceThickness,'k');
        plot3(X1*INFO.PixelSpacing(1),Y1*INFO.PixelSpacing(2),Z2*INFO.SliceThickness,'k');
    end
end
%   daspect([1,1,.3]);axis tight;
%   OptionZ.FrameRate=31;OptionZ.Duration=11;OptionZ.Periodic=true;
%   CaptureFigVid([180,90;290,45;360,0; 450, 45; 540, 90], 'WellMadeVid',OptionZ)
%     view(180,90);
end

