function [data_subpl,fitresult,time2run,outlierProps,kTPredictor,MaskVol,PathLength] = PipelineAlgorithm(fileloc,Choices,INFO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This MATLAB algorithm should comprise 13 files which correspond to Allenby's 2020 prepint on unruptured intracranial aneurysm detection, primarily from TOF MRA images.
%This code is run in MATLAB's 2018a package and is reliant on bwdistsc.m (Yuriy Mishchenko, 2007) 
%and MHAload.m (adapted Dirk Jan Kroon, 2011) functions, and several functions within MATLAB's Image Processing Toolkit. 
%please contact Mark for questions/interest. Mark.Allenby@qut.edu.au or mcallenby@gmail.com (2020 homepage: https://staff.qut.edu.au/staff/mark.allenby)

%% For examples of how to run please see 'FileLoaderExample.m'. For typical hospital scans, the three inputs are:
%'fileloc' input should be a text file of the path location to the external '.dicom' containing  folder
%'Choices' is the number corresponding to the 3D image of interest in the dataset. Find out the number using T within '[T,NewImage] = DICOMreader(fileloc)'
%'INFO' is blank e.g. 'INFO = []'. 
%These three inputs are different for other image file types ('.mha','.nii'), where the file location is the image file, choices is blank, and some INFO must be defined ( 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Collecting preliminary info.
ChoiceDisplay = 1; %if you load multiple images (in beta testing), choose which one to output.
colorstring = {'k.','g.','b.','m.','r.'}; %if you load multiple images (in beta testing), choose colors of overlapping images. Needs same length as Choices vector.

if length(Choices)==0 
    if min(fileloc((end-3):end) == '.nii')==1 
        Images{1} = niftiread(fileloc);
    elseif min(fileloc((end-3):end) == '.mha')==1 
        Images{1} = MHAload(fileloc);
    end
elseif length(INFO)==0 
    [T,NewImage] = DICOMreader(fileloc);
    for i=1:length(Choices)
        if i == 1
            [Images{i},INFO] = DICOMload(NewImage,Choices(i),fileloc);
        else
            [Images{i},~] = DICOMload(NewImage,Choices(i),fileloc);
        end
        if sum(ismember('XA',INFO.Modality))>0 %our CT images needed to be flipped and have one plane of noise deleted
            Images{i}(:,:,1) = [];
            INFO.SliceThickness = 1; 
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%All parameters below are consistent with Table 1 in Allenby et al (2020)
%manuscript, although names may have been changed for clarity. 

%%Parameter Choices
%General and Plotting Parameters
NumTissue = 2 + max(ismember('XA',INFO.Modality)); %2 for MRA, 3 for other
WhichTissue = 2 + max(ismember('XA',INFO.Modality)); %second tissue brightness for MRA, third brightness for other
ScanSpace = 7; %this dictates the size of the bounding box on the final plot, it has no other function than aesthetic.

%Threshold params
binnum = 50; %voxel intensity histogram bin number 
PeriphPrc = 0.04; %peripheral histogram bins to delete (as a percent, eg 0.04=4% is peripheral 2 bins for 50 total [3 48])
kT = 0.3134*mean(INFO.PixelSpacing')^-1.522;%This is the critical kT weight parameter which guides the extent of global or local segmentation (ratio)
    %We've assessed 3 kT weight parameters in our paper:
    %1.78 for constant relationship
    %0.6709*mean(INFO.PixelSpacing')^-0.839 for 1.20 intercept at 0.51
    %0.3134*mean(INFO.PixelSpacing')^-1.522 for 0.9 intercept at 0.51
    %0.1296*mean(INFO.PixelSpacing')^-2.222 for 0.6 intercept at 0.51

%Global reconstruction/Skeletonisation params
ObjectImageFrac = 1.3*10^-4; %the min segment object size as a fraction of total image volume, anything smaller is deleted (ratio)
SmallLinks = 0.7+max([INFO.PixelSpacing' INFO.SliceThickness]); %max distance between centrepoints of the same branch (mm)
BigLinks = 2+max([INFO.PixelSpacing' INFO.SliceThickness]); %max distance between branch connections (mm)

%%Outlier Detection
ConInt = 0.96; %the regression confidence interval beyond which outliers are considered as UIA candidates. 96% for validaiton.
minClustRate = [17.3, 1.5]/prod([INFO.PixelSpacing' INFO.SliceThickness]); %Determines minimum outlier cluster size to be a UIA candidate (mm3).
    %[17.3, 1.5] for validation; [0, 0.83 OR 0.69] for TAB 42 or other TAB 56 respectively

%%Optional Local Segmentation // Mask Refinement
reThreshIters = 1; %Rethreshold iterations. Integers of 0 and above. Might stall out with more than 2 or higher.
TightenSpacing = 0.5; %Ratio [0 1] dictates distance to travel before next resegmentation neighborhood
RegionRad = 2.75; %Resegmentation neighbourhood radius (mm)

%%Optional Mask Refinement 
RefineIA = 1; %Switch to use mask refinement (1 = yes, 0 = no)
ConIntRefine = 0.6; %The regression confidence interval beyond which outliers are considered for branch culling. [0 1] as %.
RRadMult = 10; %distance from a UIA centroid to identify a branch (mm) 
Dclose = 6.6; %distance centrepoints from end of identified branch point (mm) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Golbal threshold and align sequential images
ImageIndex = 1:length(Images);
for i = ImageIndex
    ImageHeight(i) = size(Images{i},3);
end
[~,tallImage] = max(ImageHeight);
ImageIndex(ImageIndex==tallImage) = [];
for i=1:length(Images)
    [Threshold{i},MMThresh,kTPredictor] = AUTOthresh(NumTissue,WhichTissue,kT,Images{i},binnum,ObjectImageFrac,PeriphPrc,0,1); %shuld be more like 3400 instead of 4400
    [blanki{i},data_pl{i}] = DICOMthresh(Images{i},Threshold{i},ObjectImageFrac,0);
end
MaskVol = sum(sum(sum(blanki{1})))*prod([INFO.PixelSpacing' INFO.SliceThickness]);
[optimizer, metric] = imregconfig('monomodal');
for i=ImageIndex
    tform = imregtform(blanki{i},blanki{tallImage},'affine',optimizer, metric);
    Images{i} = imwarp(Images{i},tform,'OutputView',imref3d(size(Images{i})));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For each patient checkup chosen
for i=1:length(Images)
    %% Define vascular tree and calculate geometric transforms
    blankh{i} = double(Images{i}).*double(blanki{i});
    blankd{i} = bwdistsc(imcomplement(blanki{i}),[INFO.PixelSpacing(1) INFO.PixelSpacing(2) INFO.SliceThickness]); %https://au.mathworks.com/matlabcentral/fileexchange/15455-3d-euclidean-distance-transform-for-variable-data-aspect-ratio
    RegionRad = max(max(max(blankd{1})))*RegionRad;
    [data_pl{i},blankp{i},blanks{i},blanki{i},Matched{i},BottomBranch] = FormBranchedTree(Images{i},blankd{i},blanki{i},INFO,NumTissue,WhichTissue,kT,...
        BigLinks,SmallLinks,binnum,ObjectImageFrac,MMThresh,RegionRad,TightenSpacing,PeriphPrc,reThreshIters);
    if reThreshIters>0
        blankh{i} = double(Images{i}).*double(blanki{i});
        blankd{i} = bwdistsc(imcomplement(blanki{i}),[INFO.PixelSpacing(1) INFO.PixelSpacing(2) INFO.SliceThickness]);
    end
    blanksr = bwconncomp(blanks{i},6);
    data_pls = regionprops(blanksr, 'PixelList');
    data_plsi{i} = vertcat(data_pls.PixelList);
    blankdiam{i} = bwdistsc(blanks{i},[INFO.PixelSpacing(1) INFO.PixelSpacing(2) INFO.SliceThickness]).*blanki{i};
    
    %% Define regression and identify spatially-clustered outlier regions
    [data_subpl,fitresult,blanks{i},blankp{i},blankd{i},blankdiam{i},outlierProps] = AnomalyDetection(Images{i},data_pl{i},blanks{i},blanki{i},...
        blankdiam{i},blankh{i},blankd{i},blankp{i},Matched{i},ConInt,minClustRate,SmallLinks,INFO,BottomBranch,RefineIA,ConIntRefine,Dclose,RRadMult);
    subRegions = [];
    for j = 1:length(data_subpl)
        subRegions{i,j} = [(min(data_subpl(j).PixelList(:,1))-ScanSpace),(max(data_subpl(j).PixelList(:,1))+ScanSpace),...
            (min(data_subpl(j).PixelList(:,2))-ScanSpace),(max(data_subpl(j).PixelList(:,2))+ScanSpace),...
            (min(data_subpl(j).PixelList(:,3))-ScanSpace),(max(data_subpl(j).PixelList(:,3))+ScanSpace)];
    end
    DiagPix{i} = vertcat(data_subpl.PixelList);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Metrics and figure output
MultiRegionMesh(subRegions,Images,Threshold,INFO,data_pl,blankp,blankd,blankdiam,blankh,DiagPix,data_plsi,colorstring,ChoiceDisplay, fitresult)
time2run = toc;
fitresult = [fitresult.p00,fitresult.p10,fitresult.p01,fitresult.p20,fitresult.p11,fitresult.p30,fitresult.p21,fitresult.p40,fitresult.p31];
PathLength = max(max(max(blankp{1})));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end