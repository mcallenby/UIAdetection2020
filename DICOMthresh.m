function [blanki,data_pl] = DICOMthresh(Images,Threshold,ObjectImageFrac,Voxels_min)
ImagesTEMP = Images;
if length(Threshold) == 2
    THRESH = (ImagesTEMP>=Threshold(1)).*(ImagesTEMP<=Threshold(2));
elseif length(Threshold) == 1
    THRESH = ImagesTEMP>=Threshold;
end
data_c = bwconncomp(THRESH,26); %6, 18, 26
data_pl = regionprops(data_c, 'PixelList');%, 'Area', 'BoundingBox','FilledArea');
Voxels_sub = cell2mat(cellfun(@length,cellfun(@transpose,{data_pl.PixelList},'UniformOutput',false),'uni',false));
if Voxels_min>0
    Voxels_min = cell2mat(cellfun(@(x) min(x(3,:)),cellfun(@transpose,{data_pl.PixelList},'UniformOutput',false),'uni',false)); %see if they are hooking into the ground
    index = (Voxels_sub>numel(ImagesTEMP)*ObjectImageFrac).*(Voxels_min<(size(ImagesTEMP,3)*0.2)).*(1:length(Voxels_sub)); %% 0.0001 frac is flexible for MR/CT
else
    index = (Voxels_sub>sum(sum(sum(~isnan(ImagesTEMP))))*ObjectImageFrac).*(1:length(Voxels_sub));
end
index(index == 0) = [];
blanki = zeros(size(THRESH));

if ObjectImageFrac>0
    BoundsList = zeros(1,6);
    for i = 1:length(index)
        data_pli = vertcat(data_pl(index(i)).PixelList);
        BoundsList(i,:) = [ max(data_pli,[],1)  min(data_pli,[],1) ];
    end
    if sum(BoundsList)~=0
        NoiseThreshold = 80;
        Bound2DArea = sqrt((BoundsList(:,4) - BoundsList(:,1)).^2+(BoundsList(:,5) - BoundsList(:,2)).^2);
        index(Bound2DArea'<NoiseThreshold) = [];
    end
end

data_pl = vertcat(data_pl(index).PixelList);
if ~isempty(data_pl)
    blanki(sub2ind(size(THRESH),data_pl(:,2),data_pl(:,1),data_pl(:,3))) = 1;
end

end