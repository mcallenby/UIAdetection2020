%Suggest to read comments in 'PipelineAlgorithm.m' first. This is a guide of personal example input variables for that function.

%We cannot share the RBWH DICOM dataset publically due to ethical governance, 
%but we have provided 10 example healthy MRAs from the publically avaliable 
%(but low resolution) IXI and MIDAS databases. Please go to their websites 
%for more images and information. 
%IXI: (https://brain-development.org/)
%MIDAS (https://www.insight-journal.org/midas/)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SECTION TO INDIVIDUALLY PROCESS ONE IMAGE (Load a .nii or .mha image) - uncomment below if in use.
fileloc = strcat(pwd,'\MIDAS002-MRA.mha');
DICOMnum = [];                      %Set DICOMnum to null if not a DICOM stack
INFO.PixelSpacing = [0.51; 0.51] ;  %[0.51; 0.51] MIDAS; [0.47; 0.47] IMPERIAL
INFO.SliceThickness = 0.8;          %0.8 MIDAS; 0.8 IMPERIAL
INFO.ContentDate = '';              %2010 MIDAS;  2007 (IMPERIAL);
INFO.StudyDescription = '';         %MIDAS Healthy; Imperial Healthy
INFO.Modality = 'MR';               %MR; DSA; CTA
[data_subpl,fitresult,time2run,outlierProps,KwPredictor,MaskVol,PathLength] = PipelineAlgorithm(fileloc,DICOMnum,INFO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SECTION TO INDIVIDUALLY PROCESS ONE IMAGE (Load a full .DICOM stack) - uncomment below if in use.
%fileloc = 'C:\Users\allenbym\OneDrive - Queensland University of Technology\Desktop\Aneurysm_Patients\TAB_OP_2018\DICOM\'; % (needs to be the folder containing all DICOM stacks)
%DICOMnum = [65];       %Which image number in the DICOM is desired. Can check using 'T' table output to DICOMreader (check PipelineAlgorithm).
%INFO = [];             %Set INFO to null if a DICOM stack
%[data_subpl,fitresult,time2run,outlierProps,KwPredictor,MaskVol,PathLength] = PipelineAlgorithm(fileloc,DICOMnum,INFO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SECTION TO BATCH PROCESS A FOLDER REPOSITORY (Extracted single images formatted as .MHA or .NII from a .DICOM stack) - uncomment below if in use.
% fileloc = 'U:\Research\Projects\ihbi\btm\btm_general\BTM Folder Structure\Mark_Allenby\06_ClinicalData\IntracranialAneurysms\IMPERIAL_healthyrepository\UnZipped\';
% fileloc = 'U:\Research\Projects\ihbi\btm\btm_general\BTM Folder Structure\Mark_Allenby\06_ClinicalData\IntracranialAneurysms\MIDAS_healthyrepository\';
% % % % % % % % parameters below
% PoorImages = [];
% INFO.PixelSpacing = [0.47; 0.47] ; %[0.51; 0.51] MIDAS;  [0.62; 0.62] NITRC;  [0.47; 0.47] IMPERIAL
% INFO.SliceThickness = 0.8; %0.8 MIDAS and IMPERIAL; 0.62 NITRC
% INFO.ContentDate = '2007'; %2010 MIDAS;  (NITRC);  2007 (IMPERIAL);
% INFO.StudyDescription = 'IXI Healthy'; %MIDAS Healthy; NITRC Healthy; Imperial Healthy
% INFO.Modality = 'MR';
% % % % % % % % batch loop with reporting below
% i = 1;
% j = 0;
% k = 0;
% DICOMnum = [];
% outlierProps = cell(1,1);
% data_subpl = cell(1,1);
% fitresult = [];
% time2run = [];
% Directory = dir(fileloc);
% while i<=length(Directory)
%     if sum((i+k+j) == PoorImages)>0
%         strcat('#',num2str(i+j+k),'. REMOVED. ', num2str(j), ' have been deleted so far.')
%         k = k+1;
%         RemoveList(k) = i+j+k-1;
%         Directory(i) = [];
%     elseif    length(Directory(i).name)>5
%         try
%             [data_subpl{i},fitresult(i,:),time2run(i),outlierProps{i,:},KwPredictor(i),MaskVol(i),PathLength(i)] = PipelineAlgorithm(Directory(i).name,DICOMnum,INFO);
%             UIAcounter = [];
%             UIAcounter = cell2mat(cellfun(@length,data_subpl,'uni',false));
%             Acc = sum(UIAcounter==0)/length(UIAcounter);
%             strcat('#',num2str(i+j+k),'.  Image ', num2str(i), '. UIAs: ',num2str(length(data_subpl{i})), '. Acc: ', num2str(round(Acc*100)),...
%                 '%. Kw: ',num2str(round(KwPredictor(i)*100)/100), '. MaskVol: ',num2str(round(MaskVol(i))), 'mm3. PathLength: ',num2str(round(PathLength(i))), 'mm. Time: ', num2str(round(time2run(i))), 's.')
%             i = i+1;
%         catch %this is an image with an unfixable error (2% of images have this)
%             strcat('#',num2str(i+j+k),'. DELETED. ', num2str(j+1), ' have been deleted so far.')
%             j = j+1;
%             DeleteList(j) = i+j+k-1;
%             Directory(i) = [];
%         end
%     else %this isn't an image, this is metadata text (every repository has a couple)
%         Directory(i) = [];
%     end
% end