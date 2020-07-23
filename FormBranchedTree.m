function [data_full,blankp,blanks,blanki,Matched,BottomBranch] = FormBranchedTree(Images,blankd,blanki,INFO,NumTissue,WhichTissue,Kw,BigLinks,SmallLinks,binnum,ObjectImageFrac,MMThresh,RegionRad,TightenSpacing,PeriphPrc,reThreshIters)
clear SubInd Matched

blanks = bwskel(imclose(imdilate(bwskel(imclose(imfill(blanki>0,4,'holes'),strel('sphere',1))),strel('sphere',1)),strel('sphere',1)));
for Iter = 1:reThreshIters
    [~,blanks,blanki] = CreateBranches(Images,blanks,blankd,blanki,INFO,SmallLinks,NumTissue,WhichTissue,Kw,binnum,ObjectImageFrac,MMThresh,RegionRad,PeriphPrc,TightenSpacing,1);
end
[SubInd,blanks,blanki] =CreateBranches(Images,blanks,blankd,blanki,INFO,SmallLinks,NumTissue,WhichTissue,Kw,binnum,ObjectImageFrac,MMThresh,RegionRad,PeriphPrc,TightenSpacing,0);

% ColorVec = parula(length(SubInd));
% figure
% hold on
% for k = 1:length(SubInd)
%     plot3(SubInd{k}(:,1)*INFO.PixelSpacing(1),SubInd{k}(:,2)*INFO.PixelSpacing(2),SubInd{k}(:,3)*INFO.SliceThickness,'Color',ColorVec(k,:),'LineWidth',3)
% end
% xlabel('mm')
% ylabel('mm')
% zlabel('mm')
% view(180,90);

k = 1;
while k<=length(SubInd)
    if size(SubInd{k},1)<=3
        SubInd(k) = [];
    else
        k = k+1;
    end
end

%% Puts bottom vessels in ascending, then hooks rest of vessels onto bottom in vertically sorted (bottom to top) order
ZMax = cellfun(@(x) min(x,[],1),SubInd,'UniformOutput',false);
ZMax = vertcat(ZMax{:});
[Z,ZMax] = sort(ZMax(:,3));
BottomBranch = round((max(Z)-min(Z))*0.13+min(Z));
ZInd = 1;
Matched = [];
for i = 1:size(blanki,3)%needs to be connected to the image bottom (in case of low hanging vessels) - doesnt work for LMS DSA
    blankbase = bwconncomp(blanki(:,:,i),6);
    data_base = regionprops(blankbase, 'PixelList');
    list_base = vertcat(data_base.PixelList);
    if size(list_base,1)/prod([size(Images,1) size(Images,2)])>0.0005
        break
    end
end
blanksr = bwconncomp(blanki(:,:,1:BottomBranch),6); %bwconncomp(imdilate(blanki(:,:,round((max(Z)-min(Z))*BottomFrac+min(Z))),strel('sphere',4)),6)
data_pls = regionprops(blanksr, 'PixelList');
for i = 1:length(data_pls)
    if (sum(ismember(sub2ind(size(Images(:,:,1)),data_pls(i).PixelList(:,1),data_pls(i).PixelList(:,2)),...
            sub2ind(size(Images(:,:,1)),list_base(:,1),list_base(:,2))))==0)
        data_pls(i).PixelList = [];
    end
end

[~,SInd] = sort(cell2mat(cellfun(@length,cellfun(@transpose,{data_pls.PixelList},'UniformOutput',false),'uni',false)),'descend');
i = 1;
count = zeros(size(SInd));
while i<=length(SInd)
    count(i) = 0;
    if length(data_pls(SInd(i)).PixelList)>0
        for j = 1:length(SubInd)
            if sum(ismember(sub2ind(size(Images),data_pls(SInd(i)).PixelList(:,1),data_pls(SInd(i)).PixelList(:,2),repmat(BottomBranch,...
                    size(data_pls(SInd(i)).PixelList,1),1)),sub2ind(size(Images),SubInd{j}(:,1),SubInd{j}(:,2),SubInd{j}(:,3))))>0 %0.09
                if count(i)==0
                    count(i) = j;
                elseif size(SubInd{j},1)>size(SubInd{count(i)},1)
                    count(i) = j;
                end
            end
        end
    end
    if count(i)>0
        data_pls(SInd(i)).PixelList = [];
    end
    i = i+1;
end
[~, ia, ~] = unique(count);
count = count(sort(ia,'ascend'));
count(count==0) = [];

cMatched = 0;
for i = 1:length(count)
    Matched{length(Matched)+1} = SubInd{count(i)};
    SubInd{count(i)} = [];
    Index = Matched{end}(:,3)<BottomBranch; %ismember(Matched{end}(:,3),Z(ZInd));
    if sum(Index)~=length(Index)
        ne0 = find((Index<1)~=0)';
        ix0 = unique([ne0(1) ne0(diff([0 ne0])>1)]);
        ix1 = ne0([find(diff(ne0)>1) length(ne0)]);
        if (length(ix1)==1)*(ix0(1)~=1)*(ix1(1)~=length(Index))>0 %a branch starts in one carotid and ends in another vessel at bottom (3 images of 20 have this)
            Matched{end}(ix0:ix1,4) = min([Matched{end}(ix0:ix1,4) max(Matched{end}(ix0:ix1,4))-Matched{end}(ix0:ix1,4)],[],2);
            if sum(ismember(count,count(i)))>1 %super rare case where one branch starts in one corotid and ends in another corotid (1 image of 20; TAB 42)
                cMatched = cMatched + 1;
            end
        elseif (length(ix1)==1)*(ix0(1)~=1)>0 %ascending from base %%nearly everything should be this
            Matched{end}((ix0-1):ix1,4) = Matched{end}((ix0-1):ix1,4) - Matched{end}(ix0-1,4);
        elseif (length(ix1)==1)*(ix1(1)~=length(Index))>0 %descending to end
            Matched{end}(ix0:(ix1+1),4) = Matched{end}(ix1+1,4) - Matched{end}(ix0:(ix1+1),4);
        elseif length(ix1)>1 %middle, one or multiple intersections
            for k = 1:length(ix0) %for each break
                if ix0(k) ~= 1 %it is ascending
                    Matched{end}(ix0(k):ix1(k),4) = Matched{end}(ix0(k):ix1(k),4) - min(Matched{end}((ix0(k)-1):ix1(k),4));
                elseif ix0(k)==1 %it is descending
                    Matched{end}(ix0(k):ix1(k),4) = max(Matched{end}(ix0(k):(ix1(k)+1),4))-Matched{end}(ix0(k):ix1(k),4);
                else %it is between two 0s
                    Matched{end}(ix0(k):ix1(k),4) = min([Matched{end}(ix0(k):ix1(k),4) max(Matched{end}(ix0(k):ix1(k),4))-Matched{end}(ix0(k):ix1(k),4)],[],2);
                end
            end
        end
        Matched{end}(Index,4) = Matched{end}(Index,4) - Matched{end}(Index,4);
    end
    cMatched = cMatched + 1;
%     if cMatched >= (1+max(ismember('MR',INFO.Modality))*2) %not sure whether MR should be 2 or 3... the 3rd is typically small
%         break
%     end
end

clear Index
SubRemain = cellfun(@(x) ~isempty(x),SubInd,'UniformOutput',false);
SubRemain = vertcat(SubRemain{:}).*(1:length(vertcat(SubRemain{:})))';
while sum(SubRemain)>0
    DIST = NaN(length(Matched),length(SubInd));
    Index = cell(length(Matched),length(SubInd));
    for i = 1:length(SubRemain)
        if SubRemain(i)>0
            for j = 1:length(Matched)
                Index{j,i} = ismember(sub2ind(size(Images),SubInd{i}(:,1),SubInd{i}(:,2),SubInd{i}(:,3)),...
                    sub2ind(size(Images),Matched{j}(:,1),Matched{j}(:,2),Matched{j}(:,3)));
                if sum(Index{j,i})>0
                    DIST(j,i) = min(Matched{j}(ismember(sub2ind(size(Images),Matched{j}(:,1),Matched{j}(:,2),Matched{j}(:,3)),...
                        sub2ind(size(Images),SubInd{i}(Index{j,i},1),SubInd{i}(Index{j,i},2),SubInd{i}(Index{j,i},3))),4));
                end
            end
        end
    end
    
    if sum(sum(~isnan(DIST))) == 0
        break  %%%the branches that were on base but not largest 3 are not being picked up
    end
    
    [DISTm,i] = min(DIST,[],2);
    [DIST,j] = min(DISTm,[],1);
    i = i(j);
    ne0 = find((Index{j,i}<1)~=0)';
    ix0 = unique([ne0(1) ne0(diff([0 ne0])>1)]);
    ix1 = ne0([find(diff(ne0)>1) length(ne0)]);
    
    if ((ix0(1)==1)+(ix1(1)==length(Index{j,i})))<2 %if it connects at base or if connects at tip
        Matched{length(Matched)+1} = SubInd{i};
        SubInd{i} = [];
        if ((sum(Index{j,i})==1)*(ix0(1)~=1)*(length(ix1)==1))>0 %ascending from base
            Matched{end}((ix0-1):ix1,4) = Matched{end}((ix0-1):ix1,4) - Matched{end}(ix0-1,4);
        elseif ((sum(Index{j,i})==1)*(ix1(1)~=length(Index{j,i}))*(length(ix1)==1))>0 %descending to end
            Matched{end}(ix0:(ix1+1),4) = Matched{end}(ix1+1,4) - Matched{end}(ix0:(ix1+1),4);
        elseif ((sum(Index{j,i})>1)+(length(ix1)~=1))>0 %middle, one or multiple intersections
            k=1;
            while k <= length(ix0) %for each break
                if ((ix0(k)~=1)*(ix1(k)~=length(Index{j,i})))>0 %it is between two touches
                    Matched{end}((ix0(k)-1):(ix1(k)+1),4) = min([Matched{end}((ix0(k)-1):(ix1(k)+1),4) max(Matched{end}((ix0(k)-1):(ix1(k)+1),4))-Matched{end}((ix0(k)-1):(ix1(k)+1),4)],[],2);
                elseif (ix0(k)~=1) %it is ascending
                    Matched{end}((ix0(k)-1):ix1(k),4) = Matched{end}((ix0(k)-1):ix1(k),4) - Matched{end}(ix0(k)-1,4);
                elseif (ix1(k)~=length(Index{j,i})) %it is descending
                    Matched{end}(ix0(k):(ix1(k)+1),4) = Matched{end}((ix1(k)+1),4)-Matched{end}(ix0(k):(ix1(k)+1),4);
                end
                k = k+1;
            end
        end
        Matched{end}(Index{j,i},4) = Matched{end}(Index{j,i},4) - Matched{end}(Index{j,i},4);
        Matched{end}(:,4) = Matched{end}(:,4) + DIST;
    end
    SubRemain = cellfun(@(x) ~isempty(x),SubInd,'UniformOutput',false);
    SubRemain = vertcat(SubRemain{:}).*(1:length(vertcat(SubRemain{:})))';
    clear Index DIST
end

%% Goes through all vessels and recalculates distances for complicated bifurcations - repeat makes more accurate
for counter = 1:3
    for i = 1:length(Matched)
        clear Index SUMMARY BifList
        for j = 1:length(Matched)
            Index{j} = sum(sqrt((repmat(Matched{j}(:,1),1,size(Matched{i},1)) - repmat(Matched{i}(:,1)',size(Matched{j},1),1)).^2 + ...
                (repmat(Matched{j}(:,2),1,size(Matched{i},1)) - repmat(Matched{i}(:,2)',size(Matched{j},1),1)).^2 + ...
                (repmat(Matched{j}(:,3),1,size(Matched{i},1)) - repmat(Matched{i}(:,3)',size(Matched{j},1),1)).^2) < BigLinks,2)>0;
            SUMMARY(j) = sum(Index{j});
            if i==j
                Index{j} = sqrt((Matched{j}(:,1) - Matched{i}(end,1)).^2 + ...
                    (Matched{j}(:,2) - Matched{i}(end,2)).^2 + ...
                    (Matched{j}(:,3) - Matched{i}(end,3)).^2) < BigLinks;
                SUMMARY(j) = ismember(sub2ind(size(Images),Matched{j}(end,1),Matched{j}(end,2),Matched{j}(end,3)),...
                    sub2ind(size(Images),Matched{i}(:,1),Matched{i}(:,2),Matched{i}(:,3)));
            end
        end
        BifList = (SUMMARY>0).*(1:length(SUMMARY));
        BifList(BifList==0) = [];
        if sum(SUMMARY)>0
            for k = 1:length(BifList)
                newIndex(:,k) = logical(sum(sqrt((repmat(Matched{i}(:,1),1,size(Matched{BifList(k)}(Index{BifList(k)},2),1)) - repmat(Matched{BifList(k)}(Index{BifList(k)},1)',size(Matched{i},1),1)).^2 + ...
                    (repmat(Matched{i}(:,2),1,size(Matched{BifList(k)}(Index{BifList(k)},2),1)) - repmat(Matched{BifList(k)}(Index{BifList(k)},2)',size(Matched{i},1),1)).^2 + ...
                    (repmat(Matched{i}(:,3),1,size(Matched{BifList(k)}(Index{BifList(k)},2),1)) - repmat(Matched{BifList(k)}(Index{BifList(k)},3)',size(Matched{i},1),1)).^2)<BigLinks,2)>0);
            end
            k = 1;
            while k <= length(BifList)
                if sum(newIndex(:,k))>0
                    ReMatch = Matched{i};
                    ReMatch(:,5) = NaN;
                    ne0 = find((newIndex(:,k)<1)~=0)';
                    if size(ne0,2)>0
                        ix0 = unique([ne0(1) ne0(diff([0 ne0])>1)]);
                        ix1 = ne0([find(diff(ne0)>1) length(ne0)]);
                        init = [1 ix1];
                        term = [ix0 length(newIndex(:,k))];
                        Prior = Matched{BifList(k)}(Index{BifList(k)},4) ;
                        for m = 1:length(init)
                            if ceil((m)/2)>length(Prior)
                                Prior = [Prior,Prior];
                            end
                            if init(m)~=term(m)
                                ReMatch(init(m):term(m),5) = Prior(ceil((m)/2)) + min(sqrt((repmat(Matched{i}((init(m)):term(m),1),[1,size(Matched{BifList(k)}(Index{BifList(k)},1))]) - repmat(Matched{BifList(k)}(Index{BifList(k)},1)',[size(Matched{i}((init(m)):term(m),1)),1])).^2+...
                                    (repmat(Matched{i}((init(m)):term(m),2),[1,size(Matched{BifList(k)}(Index{BifList(k)},1))]) - repmat(Matched{BifList(k)}(Index{BifList(k)},2)',[size(Matched{i}((init(m)):term(m),1)),1])).^2+...
                                    (repmat(Matched{i}((init(m)):term(m),3),[1,size(Matched{BifList(k)}(Index{BifList(k)},1))]) - repmat(Matched{BifList(k)}(Index{BifList(k)},3)',[size(Matched{i}((init(m)):term(m),1)),1])).^2),[],2);
                            end
                        end
                        if ((sum(newIndex(:,k))==1)*(ix0(1)~=1)*(length(ix1)==1))>0 %ascending from base
                            ReMatch((ix0-1):ix1,5) = min(Matched{BifList(k)}(Index{BifList(k)},4)) + abs(cumsum([0;Matched{i}((ix0(1)-1):(ix1(1)-1),4) - Matched{i}(ix0(1):ix1(1),4)]));
                        elseif ((sum(newIndex(:,k))==1)*(ix1(1)~=length(newIndex(:,k)))*(length(ix1)==1))>0 %descending to end
                            ReMatch(ix0:(ix1+1),5) = min(Matched{BifList(k)}(Index{BifList(k)},4)) + abs(flipud(cumsum(flipud([Matched{i}(ix0(1):(ix1(1)),4) - Matched{i}((ix0(1)+1):(ix1(1)+1),4); 0]))));
                        elseif ((sum(newIndex(:,k))>1)+(length(ix1)~=1))>0 %middle, one or multiple intersections
                            m=1;
                            while m <= length(ix0) %for each break
                                if ((ix0(m)~=1)*(ix1(m)~=length(newIndex(:,k))))>0 %it is between two touches
                                    ReMatch((ix0(m)-1):(ix1(m)+1),5) = min([ReMatch(ix1(m)+1,5) + flipud(cumsum(flipud([abs(Matched{i}((ix0(m)):(ix1(m)+1),4)-Matched{i}((ix0(m)-1):(ix1(m)),4));0]))), ...
                                        ReMatch(ix0(m)-1,5) + cumsum(abs([0;Matched{i}((ix0(m)):(ix1(m)+1),4)-Matched{i}((ix0(m)-1):(ix1(m)),4)]))],[],2);
                                elseif (ix0(m)~=1) %it is ascending
                                    ReMatch((ix0(m)-1):ix1(m),5) = ReMatch(ix0(m)-1,5)  + cumsum([0;abs(Matched{i}((ix0(m)-1):(ix1(m)-1),4)-Matched{i}((ix0(m)):(ix1(m)),4))]);
                                elseif (ix1(m)~=length(newIndex(:,k))) %it is descending
                                    ReMatch(ix0(m):(ix1(m)+1),5) = ReMatch(ix1(m)+1,5) + flipud(cumsum(flipud(abs([Matched{i}((ix0(m)+1):(ix1(m)+1),4)-Matched{i}((ix0(m)):(ix1(m)),4);0]))));
                                end
                                m = m+1;
                            end
                        end
                    else
                        ReMatch(:,5) = repmat(min(Matched{BifList(k)}(Index{BifList(k)},4)),size(Matched{i}(:,4),1),1);
                    end
                    Matched{i}(:,4) = min(ReMatch(:,4:5),[],2);
                end
                k = k+1;
            end
            newIndex(:,BifList == i) = [];
            BifList(BifList == i) = [];
            clear ReMatch
            clear newIndex
        end
        clear Index
    end
end

%  ColorVec = parula(length(Matched));
%  figure
%  hold on
%  for k = 1:length(Matched)
%      plot3(Matched{k}(:,1)*INFO.PixelSpacing(1),Matched{k}(:,2)*INFO.PixelSpacing(2),Matched{k}(:,3)*INFO.SliceThickness,'Color',ColorVec(k,:),'LineWidth',3)
%  end
%      xlabel('mm')
%      ylabel('mm')
%      zlabel('mm')
%      view(180,90);


%% hook to individual pixels and plot
FullList = vertcat(Matched{:});

blanksr = bwconncomp(blanki,6);
data_parts = regionprops(blanksr, 'PixelList');
i=1;
while i<=length(data_parts)
    if sum(ismember(sub2ind(size(Images),FullList(:,1),FullList(:,2),FullList(:,3)),sub2ind(size(Images),data_parts(i).PixelList(:,1),data_parts(i).PixelList(:,2),data_parts(i).PixelList(:,3))))==0
        data_parts(i) = [];
    else
        i = i+1;
    end
end
data_fullNEW = vertcat(data_parts.PixelList);
blanki = zeros(size(blanki));
blanki(sub2ind(size(Images),data_fullNEW(:,2),data_fullNEW(:,1),data_fullNEW(:,3))) = 1;

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

%figure
%scatter3(data_full2(:,1),data_full2(:,2),data_full2(:,3),1,data_full2(:,4))%,'MarkerEdgeAlpha', 0.05,'MarkerFaceAlpha', 0.05)
figure
scatter3(FullList(:,1)*INFO.PixelSpacing(1),FullList(:,2)*INFO.PixelSpacing(2),FullList(:,3)*INFO.SliceThickness,10,FullList(:,4),'o','filled')
    xlabel('mm')
    ylabel('mm')
    zlabel('mm')
    view(180,90);
% hold on
%scatter3(Matched{i}(:,1),Matched{i}(:,2),Matched{i}(:,3),250,Matched{i}(:,4),'.')
pause(0.01)
end