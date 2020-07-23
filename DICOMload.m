function [Images,INFO] = DICOMload(NewImage,Choice,fileloc)

ImageSelection = [NewImage(Choice),NewImage(Choice+1)-1];
i = 1;
filename = sprintf('I%d', ImageSelection(1));
INFO = dicominfo(strcat(fileloc,filename));
[dimX, dimY] = size(dicomread(strcat(fileloc,filename)));
Images = repmat(int16(0), [dimX dimY (ImageSelection(2) - ImageSelection(1))]);
for p=(ImageSelection(1)):ImageSelection(2)
    filename = sprintf('I%d', p);
    Images(:,:,i) = dicomread(strcat(fileloc,filename));
    i = i+1;
end

end
