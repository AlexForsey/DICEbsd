function  [ebsd,grains,DIC] = plotEBSD_DIC(ebsd,grains,DIC,filterVal,deform_switch)
%plotting
setMTEXpref('showMicronBar','off')
setMTEXpref('showCoordinates','on')

cropControl = questdlg('Auto crop?', 'Crop control', 'auto', 'manual', 'auto');

switch cropControl
    case 'auto'
        %crop data to extremes of the DIC strain data
        xcrop=[max(DIC.posX(:)),min(DIC.posX(:))];
        ycrop=[max(DIC.posY(:)),min(DIC.posY(:))];
        ebsd=ebsd(ebsd.x>xcrop(1) & ebsd.x<xcrop(2) & ebsd.y>ycrop(1) & ebsd.y<ycrop(2));
        [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',5*degree);
        
        h = figure(1);
        plot(ebsd,ebsd.Exx)
        hold on
        plot(grains.boundary,'linewidth',2)
        caxis([0 0.3])
        mtexColorMap hot
        hold off
        title('Exx')
        set(h,'position',[377 17 710 533])
        
    case 'manual'
        %manual crop of data
        cropLoop = 'no';
        while strcmp(cropLoop,'no')
            %figure to plot from
            h = figure(1);
            plot(ebsd,ebsd.Exx)
            caxis([0 0.3])
            hold on
            plot(grains.boundary,'linewidth',2)
            mtexColorMap hot
            hold off
            title('Click opposing corners to crop plots')
            set(h,'position',[377 17 710 533])
            
            %crop data for plotting
            [xcrop,ycrop]= ginput(1);
            hold on
            plot(xcrop,ycrop,'or','LineWidth',2)
            [xcrop(2),ycrop(2)]= ginput(1);
            plot(xcrop(2),ycrop(2),'or','LineWidth',2)
            hold off
            xcrop=sort(xcrop);y=sort(ycrop);
            %     ebsd=ebsd(ebsd.x>x(1) & ebsd.x<x(2) & ebsd.y>y(1) & ebsd.y<y(2));
            ebsd2=ebsd(ebsd.x>xcrop(1) & ebsd.x<xcrop(2) & ebsd.y>ycrop(1) & ebsd.y<ycrop(2));
            [grains2,ebsd2.grainId,ebsd2.mis2mean] = calcGrains(ebsd2,'angle',5*degree);
            
            plot(ebsd2,ebsd2.Exx)
            %     plot(ebsd,ebsd.Exx)
            hold on
            plot(grains2.boundary,'linewidth',2)
            caxis([0 0.3])
            mtexColorMap hot
            hold off
            title('Exx')
            
            cropLoop = questdlg('Accept crop?', 'Crop control', 'yes', 'no', 'yes');
        end
        grains = grains2;
        ebsd = ebsd2;
        clear grains2 ebsd2

end

if isfield(DIC,'raw')
    [grains,DIC] = meanGrainStrain(ebsd,grains,DIC,deform_switch,filterVal);
else
    grains = grainStrain(ebsd,grains);
end

h = figure(2);
plot(ebsd,ebsd.Exy)
hold on
plot(grains.boundary,'LineWidth',2)
caxis([-0.2 0.2])
mtexColorMap jet
hold off
title('Exy')
set(h,'position',[377 17 710 533])

h = figure(3);
plot(ebsd,ebsd.Eyx)
hold on
plot(grains.boundary,'LineWidth',2)
caxis([-0.2 0.2])
mtexColorMap jet
hold off
title('Eyx')
set(h,'position',[377 17 710 533])

h = figure(4);
plot(ebsd,ebsd.Eyy)
hold on
plot(grains.boundary,'LineWidth',2)
caxis([-.1 0.05])
mtexColorMap jet
hold off
title('Eyy')
set(h,'position',[377 17 710 533])

h = figure(5);
plot(ebsd,ebsd.EpMax)
hold on
plot(grains.boundary,'LineWidth',2)
%caxis([-.1 0.05])
mtexColorMap jet
hold off
title('Emax')
set(h,'position',[377 17 710 533])

h = figure(6);
plot(ebsd,ebsd.EpMin)
hold on
plot(grains.boundary,'LineWidth',2)
%caxis([-.1 0.05])
mtexColorMap jet
hold off
title('Emin')
set(h,'position',[377 17 710 533])

h = figure(7);
plot(ebsd,ebsd.Gmax)
hold on
plot(grains.boundary,'LineWidth',2)
%caxis([-.1 0.05])
mtexColorMap jet
hold off
title('Shear max')
set(h,'position',[377 17 710 533])

h = figure(8);
plot(ebsd,ebsd.Theta)
hold on
plot(grains.boundary,'LineWidth',2)
%caxis([-.1 0.05])
mtexColorMap jet
hold off
title('Principal strain orientation')
set(h,'position',[377 17 710 533])

h = figure(9);
plot(ebsd,ebsd.Rotation)
hold on
plot(grains.boundary,'LineWidth',2)
%caxis([-.1 0.05])
mtexColorMap jet
hold off
title('rigid body rotation angle')
set(h,'position',[377 17 710 533])

h = figure(10);
plot(grains,grains.Exx)
hold on

plot(grains.boundary,'LineWidth',2)
hold off
mtexColorMap hot
caxis([0 0.2])
title('Average Exx in grains')
set(h,'position',[377 17 710 533])

h = figure(11);
plot(ebsd,ebsd.Exx)
hold on
plot(grains.boundary,'LineWidth',2)
mtexColorMap hot
hold off
caxis([0 0.2])
title('Average Exx in grains')
set(h,'position',[377 17 710 533])

end