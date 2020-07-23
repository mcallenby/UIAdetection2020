function [ThresholdFinal,MMThresh,KwWeight] = AUTOthresh(NumTissue,WhichTissue,KwWeight,Images,binnum,ObjectImageFrac,PeriphPrc,MMThresh,PLOTON)
Images(Images<0) = 0;
[y,edges] = histcounts(double(Images(:)),binnum);
x = (edges(2:end) - edges(1:(end-1)))/2+ edges(1:(end-1));
y(y==0) = 1;

fvalstore = 1000;
pinf = [];
iter2 = 0;

BaseAcc = 0.04;
IterAccDecay = 0.005; 
IterNum = 25;
while (fvalstore>(BaseAcc+IterAccDecay*iter2)) || (length(pinf)<NumTissue)
    iter2 = iter2+1;
    p0 = [logspace(log10(sum(log(y))*10^2),log10(sum(log(y))*1.5*10^1),NumTissue) logspace(log10(1.5),log10(0.5),NumTissue) logspace(log10(edges(length(edges))/3),log10(edges(end-4)),NumTissue)];
    A = [];
    b = []; 
    Aeq = [];
    beq = [];
    lb = [logspace(log10(sum(log(y))*10^1),log10(sum(log(y))*10^0),NumTissue) logspace(log10(0.1),log10(0.01),NumTissue) logspace(log10(edges(2)),log10(edges(length(edges))/2),NumTissue)];
    ub = [logspace(log10(sum(log(y))*10^2),log10(sum(log(y))*1.5*10^1),NumTissue) logspace(log10(1.5),log10(0.5),NumTissue) logspace(log10(edges(length(edges))/3),log10(edges(end-4)),NumTissue)];
     options = optimoptions('fmincon','Display','off');
    iter = 0;
    while (iter < IterNum) && (fvalstore>BaseAcc)
        iter = iter+1;
        [p,fval,~,~] = fmincon(@(p)LogNormal(y,x,p,NumTissue,round(binnum*PeriphPrc)),p0,A,b,Aeq,beq,lb,ub,[],options);
        if fval<fvalstore
            fvalstore = fval;
            pstore = p;
        end
        p0 = [p0(1:NumTissue).*abs(randn(1,NumTissue)).^2, p0((NumTissue*2-(NumTissue-1)):(NumTissue*2)).*rand(1,NumTissue),...
            p0((NumTissue*3-(NumTissue-1)):(NumTissue*3)).*rand(1,NumTissue)];
    end
    c = pstore(1:NumTissue)';
    sig = pstore((NumTissue*2-(NumTissue-1)):(NumTissue*2))';
    mode = pstore((NumTissue*3-(NumTissue-1)):(NumTissue*3))';
    mu = log(mode)+sig.^2;
    xplot = edges(1):((edges(end)-edges(1))/1000):edges(end);
    PD = c./(sqrt(2*pi).*repmat(xplot,NumTissue,1).*sig) .* exp(-((log(repmat(xplot,NumTissue,1))-mu).^2)./(2*sig.^2));
    PDd2 = diff(sum(PD,1),2);
    ne0 = find((PDd2<0)~=0);
    pinf = [ne0(find(diff(ne0)>1)) ne0(end)+1];
    ne1 = find((PDd2>0)~=0);
    ninf = ne1(find(diff(ne1)>1));
    PDd3 = diff(sum(PDd2,1),1);
    ne1 = find((PDd3>0)~=0);
    ninf3 = ne1(find(diff(ne1)>1));
    if iter2>10
        break
    end
end

for i = 1:NumTissue
    [~, whichpeak(i)] = min(abs(xplot - mode(i)));
end

Threshold = [];
for i = 1:(NumTissue+1)
    if i==1
        if ~isempty(ninf)==1
            if min(ninf)<min(pinf)
                Threshold(i) = max([xplot(1) xplot(min(ninf))]);
            end
        end
        if isempty(Threshold)
            Threshold(i) = xplot(1);
        end
    elseif i>length(pinf)
        Threshold(i) = min([xplot(pinf(end)) xplot(end)]);
    else
        Threshold(i) = xplot(min(ninf3(xplot(ninf3)>mode(i-1))))+(mode(i)-xplot(min(ninf3(xplot(ninf3)>mode(i-1)))))*(PD(2,whichpeak(2))/(PD(1,whichpeak(1))+PD(2,whichpeak(2)))).^KwWeight; %1.2, 3.7 before
    end
end

for i = 1:(NumTissue+1)
    if i==WhichTissue
        ThresholdFinal = Threshold(i);
    elseif i==(WhichTissue+1)
        if i<(NumTissue+1)
            ThresholdFinal(2) = Threshold(i);
        end
    end
end
if MMThresh==0
    if WhichTissue==1
        MMThresh(1) = xplot(1);
    else
        MMThresh(1) = pinf(WhichTissue-1);
    end
    if WhichTissue==NumTissue
        MMThresh(2) = xplot(end);
    else
        MMThresh(2) = mode(WhichTissue+1); %better this or better pinf(WhichTissue+1)?
    end
end
ThresholdFinal(ThresholdFinal<MMThresh(1)) = MMThresh(1); %needs to be universal from first
ThresholdFinal(ThresholdFinal>MMThresh(2)) = MMThresh(2);

if PLOTON == 1
    figure
    h = histogram(double(Images(:)),binnum);
    set(gca,'YScale','log')
    hold on
    plot(xplot,exp(PD(1,:)),'Color',[0.8 0 0.8],'LineWidth',2)
    plot(xplot,exp(PD(2,:)),'Color',[0.91 0.41 0.17],'LineWidth',2)
    if NumTissue==3
        plot(xplot,exp(PD(2,:)),'Color',[0.91 0.41 0.17],'LineWidth',2)
        plot(xplot,exp(PD(3,:)),'Color',[0.8 0.8 0],'LineWidth',2)
    else
        plot(xplot,exp(PD(2,:)),'Color',[0.8 0.8 0],'LineWidth',2)
    end
    plot(xplot,exp(sum(PD,1)),'k','LineWidth',2)
    
    ThresholdFinal = [Threshold(1) Threshold(2)];
    YLims = get(gca,'YLim');
    patch('vertices', [ThresholdFinal(1), YLims(1); ThresholdFinal(1), YLims(2); ThresholdFinal(2), YLims(2); ThresholdFinal(2) YLims(1)],...
        'faces', [1, 2, 3, 4],'FaceColor', [ 0.8 0 0.8],'FaceAlpha', 0.2,'LineWidth',2,'LineStyle','--','EdgeColor',[ 0.8 0 0.8])
    
    if NumTissue==3
        ThresholdFinal = Threshold(2:3);
        patch('vertices', [ThresholdFinal(1), YLims(1); ThresholdFinal(1), YLims(2); ThresholdFinal(2), YLims(2); ThresholdFinal(2) YLims(1)],...
            'faces', [1, 2, 3, 4],'FaceColor', [ 0.91 0.41 0.17],'FaceAlpha', 0.2,'LineWidth',2,'LineStyle','--','EdgeColor',[ 0.91 0.41 0.17])
        ThresholdFinal = [Threshold(3) xplot(end)];
        patch('vertices', [ThresholdFinal(1), YLims(1); ThresholdFinal(1), YLims(2); ThresholdFinal(2), YLims(2); ThresholdFinal(2) YLims(1)],...
            'faces', [1, 2, 3, 4],'FaceColor', [ 0.8 0.8 0],'FaceAlpha', 0.2,'LineWidth',2,'LineStyle','--','EdgeColor',[ 0.8 0.8 0])
    else
        ThresholdFinal = [Threshold(2) xplot(end)];
        patch('vertices', [ThresholdFinal(1), YLims(1); ThresholdFinal(1), YLims(2); ThresholdFinal(2), YLims(2); ThresholdFinal(2) YLims(1)],...
            'faces', [1, 2, 3, 4],'FaceColor', [ 0.8 0.8 0],'FaceAlpha', 0.2,'LineWidth',2,'LineStyle','--','EdgeColor',[ 0.8 0.8 0])
    end
    
    
    [blank1,~] = DICOMthresh(Images,[xplot(ninf(1)) Threshold(2)],ObjectImageFrac,0);
    [blank2,~] = DICOMthresh(Images,[Threshold(2) xplot(end)],ObjectImageFrac,0);
    if NumTissue==3
        [blank3,~] = DICOMthresh(Images,Threshold(3:4),ObjectImageFrac,0);
    end
    
    orange = cat(3, ones(size(Images(:,:,1)))*0.91,ones(size(Images(:,:,1)))*0.41, ones(size(Images(:,:,1)))*0.17);
    yellow = cat(3, ones(size(Images(:,:,1)))*0.8,ones(size(Images(:,:,1)))*0.8, zeros(size(Images(:,:,1))));
    purple = cat(3, ones(size(Images(:,:,1)))*0.8,zeros(size(Images(:,:,1))), ones(size(Images(:,:,1)))*0.8);
    
    figure
    imshow(mean((Images),3)/double(max(max(max(Images))))*2)
    hold on
    f = imshow(purple);
    set(f,'AlphaData', mean((blank1),3)*0.1)
    if NumTissue==3
        h = imshow(orange);
        set(h,'AlphaData', mean((blank2),3))
        f = imshow(yellow);
        set(f,'AlphaData', mean((blank3),3)*10)
    else
        f = imshow(yellow);
        set(f,'AlphaData', mean((blank2),3)*10)
    end
end
pause(1)
end

function [error,x0] = LogNormal(y,x,x0,NumTissue,PeriphBins)

c = x0(1:NumTissue)';
sig = x0((NumTissue*2-(NumTissue-1)):(NumTissue*2))';
mode = x0((NumTissue*3-(NumTissue-1)):(NumTissue*3))'; %this is the mode

mu = log(mode)+sig.^2; %this is the actual mu in the equation

PD = c./(sqrt(2*pi).*repmat(x,NumTissue,1).*sig) .* exp(-((log(repmat(x,NumTissue,1))-mu).^2)./(2*sig.^2));
PD = sum(PD,1);

error = sum(abs(PD((1+PeriphBins):(end-PeriphBins))-log(y((1+PeriphBins):(end-PeriphBins))))/log(y((1+PeriphBins):(end-PeriphBins))));

end
