function DICEplot3(handles)
%define bins for hist plots
grainHist = 20;
ebsdHist = 100;
Cax(1) = str2double(get(handles.caxmin,'string'));
Cax(2) = str2double(get(handles.caxmax,'string'));

colcs=get(handles.colorMap_m_pop,'string');choi=get(handles.colorMap_m_pop,'value');
%colcs=get(handles.colorMap_m,'string');choi=get(handles.colorMap_m,'value');
colc=colcs{choi};if choi==6; colc='LaboTeX';end

plotPosSwitch = get(handles.docPlot,'Value');
switch plotPosSwitch
    case 0
        figure;h = gca;
        
        posa(1)=str2double(get(handles.p1_t,'string'));posa(2)=str2double(get(handles.p2_t,'string'));
        posa(3)=str2double(get(handles.p3_t,'string'));posa(4)=str2double(get(handles.p4_t,'string'));
    case 1
        h = handles.axes1;
        cla(h)
%         set(h,'NextPlot','replacechildren')
        setMTEXpref('xAxisDirection','west');
        setMTEXpref('zAxisDirection','inToPlane');
end

%load data to plot
plotting = getappdata(handles.output,'plotting');
load([plotting.ProjFilePathName,plotting.ProjFileName])

%select what to plot
plotVal = get(handles.plotType,'Value');
plotStr = get(handles.plotType,'string');
plotSwitch = cell2mat(plotStr(plotVal));
%select what to plot(x) for scatter
xAxisVal = get(handles.xAxisVar,'Value');
xAxisStr = get(handles.xAxisVar,'string');
xAxisSwitch = cell2mat(xAxisStr(xAxisVal));

%select type of plot
modVal = get(handles.plotMod,'Value');
modStr = get(handles.plotMod,'string');
modSwitch = cell2mat(modStr(modVal));

%Define item to plot

%make displayed propety match the variable name
if strcmp(plotSwitch,'Grain Size')
    plotSwitchField = grainSize;
elseif strcmp(plotSwitch,'Grain ID')
    plotSwitchField = 'id';
elseif strcmp(plotSwitch,'Max Principal Strain')
    plotSwitchField = 'EpMax';
elseif strcmp(plotSwitch,'Min Principal Strain')
    plotSwitchField = 'EpMin';
elseif strcmp(plotSwitch,'Angle Principal Strain')
    plotSwitchField = 'Theta';
elseif strcmp(plotSwitch,'Max Shear Strain')
    plotSwitchField = 'Gmax'; 
elseif strcmp(plotSwitch,'grainID')
    plotSwitchField = 'grainId';
else
    plotSwitchField = plotSwitch;
end

%separate value for plotting
ooo = 'o';
if isfield(ebsd,plotting)
    ebsd.plotting = [];
    eval(['ebsd.plotting = ebsd.', plotSwitchField])
else
    eval(['ebsd.prop.plotting = ebsd.', plotSwitchField])
end

if isfield(grains,plotting)
    grains.plotting = [];
    eval(['grains.plotting = grains.', plotSwitchField])
else
    eval(['grains.prop.plotting = grains.', plotSwitchField])
end


%calc extra parameters
% % [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',5*degree);

if exist('ebsd','var') && exist('grains','var')
%     switch plotSwitch
%         case 'Exx'
            switch modSwitch
                case 'Standard'
                    switch plotPosSwitch
                        case 0
                            plot(ebsd,ebsd.plotting);
                        case 1
                            plot(ebsd,ebsd.plotting,'parent',h);
                    end
                case 'Averaged'
                    switch plotPosSwitch
                        case 0
                            plot(grains,grains.plotting);
                        case 1
                            plot(grains,grains.plotting,'parent',h);
                    end
                case 'COMP: Max Schmid'
                    plotDICs(1,grains,grains,2,h);
                case 'Hist: Grains'
%                     plotDICs(1,grains,grains,1,h);
                    hist(h,grains.plotting(grains.plotting>Cax(1) & grains.plotting<Cax(2)),grainHist);
                    xlabel(plotSwitch)
                    ylabel('Frequency')
                case 'Hist: Points'
                    hist(h,ebsd.plotting(ebsd.plotting>Cax(1) & ebsd.plotting<Cax(2)),ebsdHist);
                    xlabel(plotSwitch)
                    ylabel('Frequency')
%                     hist(h,ebsd.Exx(:),100);
                case 'GOS'
                    plot(grains.prop.plotting,grains.GOS,'o')
%                     plotDICs(1,grains,grains,3,h);
                case 'Grain Size'
                    plot(grains.prop.plotting,grains.grainSize,'o')
%                     plotDICs(1,grains,grains,4,h);
                case 'Grain ID'
                    plot(grains.prop.plotting,grains.id,'o')
                case 'Scatter'
                    switch xAxisSwitch
                        case 'Exx'
                            plot(ebsd.plotting(:),ebsd.Exx(:),'o')
                            xlabel(plotSwitch);ylabel('Exx');
                        case 'Eyy'
                            plot(ebsd.plotting(:),ebsd.Eyy(:),'o')
                            xlabel(plotSwitch);ylabel('Eyy');
                        case 'Exy'
                            plot(ebsd.plotting(:),ebsd.Exy(:),'o')
                            xlabel(plotSwitch);ylabel('Exy');
                        case 'Eyx'
                            plot(ebsd.plotting(:),ebsd.Eyx(:),'o')
                            xlabel(plotSwitch);ylabel('Eyx');
                        case 'Max Principal Strain'
                            plot(ebsd.plotting(:),ebsd.EpMax(:),'o')
                            xlabel(plotSwitch);ylabel('Max principal stress');
                        case 'Min Principal Strain'
                            plot(ebsd.plotting(:),ebsd.EpMin(:),'o')
                            xlabel(plotSwitch);ylabel('Min principal stress');
                        case 'Angle Principal Strain'
                            plot(ebsd.plotting(:),ebsd.Theta(:),'o')
                            xlabel(plotSwitch);ylabel('Maximum shear angle');
                        case 'Max Shear Strain'
                            plot(ebsd.plotting(:),ebsd.Gmax(:),'o')
                            xlabel(plotSwitch);ylabel('Maximum shear strain');
                        case 'Rotation'
                            plot(ebsd.plotting(:),ebsd.Rotation(:),'o')
                            xlabel(plotSwitch);ylabel('Rotation');
                        case 'E12'
                            plot(ebsd.plotting(:),(ebsd.Exy(:) + ebsd.Eyx(:))/2,'o')
                            xlabel(plotSwitch);ylabel('E12');
                        case 'KAM'
                            plot(ebsd.plotting(:),ebsd.KAM(:),'o');
                            xlabel(plotSwitch);ylabel('KAM');
                        case 'mis2mean'
                            plot(ebsd.plotting(:),ebsd.mis2mean(:),'o');
                            xlabel(plotSwitch);ylabel('mis2mean');
                        case 'grainID'
                            plot(ebsd.plotting(:),ebsd.grainId(:),'o');
                            xlabel(plotSwitch);ylabel('grainID');
                        case 'Poisson'
                            plot(ebsd.plotting(:),ebsd.Poisson(:),'o')
                            xlabel(plotSwitch);ylabel('Poisson');
                    end
            end
            title(plotSwitch)
else
    %plot just DIC data
    if min(isfield(DIC,{'Ep' 'Gmax' 'Theta' 'RotMat'})) == 0
        for i = 1:size(DIC.E,3)
            for j = 1:size(DIC.E,4)
                [DIC.Ep(:,:,i,j),DIC.Gmax(i,j),DIC.Theta(i,j),DIC.RotMat(:,:,i,j)] = principalStrains(squeeze(DIC.E(:,:,i,j)));
            end
        end
        save([plotting.ProjFilePathName,plotting.ProjFileName],'DIC')
    end
    %     switch modSwitch
    %                 case 'Standard'
    switch plotSwitch
        case 'Exx'
            imagesc(DIC.Eposx(1,:),DIC.Eposy(:,1),squeeze(DIC.E(1,1,:,:)));
            title('Exx')
        case 'Exy'
            imagesc(DIC.Eposx(1,:),DIC.Eposy(:,1),squeeze(DIC.E(1,2,:,:)));
            title('Exy')
        case 'Eyx'
            imagesc(DIC.Eposx(1,:),DIC.Eposy(:,1),squeeze(DIC.E(2,1,:,:)));
            title('Eyx')
        case 'Eyy'
            imagesc(DIC.Eposx(1,:),DIC.Eposy(:,1),squeeze(DIC.E(2,2,:,:)));
            title('Eyx')
        case 'Max Shear Strain'
            imagesc(DIC.Eposx(1,:),DIC.Eposy(:,1),squeeze(DIC.Ep(1,1,:,:))-squeeze(DIC.Ep(2,2,:,:)));
            title('Max Shear Strain')
        case 'Max Principal Strain'
            imagesc(DIC.Eposx(1,:),DIC.Eposy(:,1),squeeze(DIC.Ep(1,1,:,:)));
            title('Max Principal Strain')
        case 'Min Principal Strain'
            imagesc(DIC.Eposx(1,:),DIC.Eposy(:,1),squeeze(DIC.Ep(2,2,:,:)));
            title('Min Principal Strain')
        case 'Angle Principal Strain'
            imagesc(DIC.Eposx(1,:),DIC.Eposy(:,1),squeeze(DIC.Theta));
            title('Angle Principal Strain')
        case 'Rotation'
            imagesc(DIC.Eposx(1,:),DIC.Eposy(:,1),squeeze(DIC.W(1,2,:,:)));
            title('Rotation')
        case 'E12'
            imagesc(DIC.Eposx(1,:),DIC.Eposy(:,1),(squeeze(DIC.E(1,2,:,:))+squeeze(DIC.E(2,1,:,:)))/2);
            title('E12')
    end
% end
    axis image
    mtexColorbar
end
%set colourmap
mtexColorMap(colc)

if exist('ebsd','var') && exist('grains','var')
    if strcmp(modSwitch,'Standard')||strcmp(modSwitch,'Averaged')
        switch plotPosSwitch
            case 0
                hold on
                plot(grains.boundary,'linewidth',2);
                hold off
            case 1
                %hold(h,'on')
                plot(grains.boundary,'parent',h,'linewidth',2);
                %hold(h,'off')
                %set(h,'XDir','reverse')
                set(h,'XDir','normal')
                axis image
        end
        mtexColorbar
        
    else
        mtexColorbar('off')
    end
else
    mtexColorbar
end
% set(h,'xdir','reverse')
% set(h,'ydir','reverse')
end

function [e1,Gmax,Theta,W] = principalStrains(e)
%Separates rigid body rotation matrix and strain tensor. Calculates min and
%max principal strains and principal strain angle.

%Alex Forsey    Apr 2016

%calculate rigid body rotation
W = (e-e')/2;

%calculate strain tensor with rigid body rotation removed
ee = e - W;

%calculate principal strain direction
Theta = (atan((ee(1,2)+ee(2,1))/2)/(ee(1,1)-ee(2,2)))/2;

%calculate transformation matrix
Q = [cos(Theta),sin(Theta);-sin(Theta),cos(Theta)];

%calculate principal strains
e1 = Q*ee*Q';

%calculate max shear strain
Gmax = max(e1(:)) - min(e1(:));

end