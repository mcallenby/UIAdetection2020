function [T,NewImage, INFO] = DICOMreader(fileloc)

DICOMnum = length(dir(fileloc));

%Creates an image index of the first image of a series
filename = sprintf('I%d', 0);
INFO = dicominfo(strcat(fileloc,filename));
j = 1;
NewImage = 0;
StackLength = [];
clear PRIORDESC
try
    PRIORDESC=cellstr(INFO.StudyDescription);
catch
    PRIORDESC = {'Blank'};
end
DETAIL=cellstr(INFO.StudyDate);
try
    PRIORDATE = cellstr(INFO.AcquisitionDate);
catch
    PRIORDATE = cellstr(INFO.StudyDate);
end
PRIORDIM = [INFO.Width,INFO.Height];
try
    for p=1:(DICOMnum-1)
        filename = sprintf('I%d', p);
        INFO = dicominfo(strcat(fileloc,filename));
        if isfield(INFO,'SeriesDescription')>0
            if (isequal(DETAIL{end},INFO.SeriesDescription)<1) || (isequal([INFO.Width,INFO.Height],PRIORDIM)<1)
                StackLength = [StackLength,p-NewImage(end)];
                NewImage = [NewImage,p];
                PRIORDESC = [PRIORDESC,cellstr(INFO.StudyDescription)];
                DETAIL = [DETAIL,INFO.SeriesDescription];
                try
                    PRIORDATE = [PRIORDATE,cellstr(INFO.AcquisitionDate)];
                catch
                    PRIORDATE = [PRIORDATE,cellstr(INFO.StudyDate)];
                end
                PRIORDIM = [INFO.Width,INFO.Height];
                j = j+1;
            end
        else
            if (isequal(PRIORDESC{end},INFO.StudyDescription)<1) || (isequal([INFO.Width,INFO.Height],PRIORDIM)<1)
                StackLength = [StackLength,p-NewImage(end)];
                NewImage = [NewImage,p];
                PRIORDESC = [PRIORDESC,cellstr(INFO.StudyDescription)];
                DETAIL = [DETAIL,'NULL'];
                try
                    PRIORDATE = [PRIORDATE,cellstr(INFO.AcquisitionDate)];
                catch
                    PRIORDATE = [PRIORDATE,cellstr(INFO.StudyDate)];
                end
                PRIORDIM = [INFO.Width,INFO.Height];
                j = j+1;
            end
        end
    end
end
StackLength = [StackLength,p-NewImage(end)];


T = table([1:j]',PRIORDATE',PRIORDESC',DETAIL',StackLength');
T.Properties.VariableNames = {'Num','Date','Type','Desc','Stacks'};

end
