function [ebsd,grains,DIC] = straininterp7_7GUI(ProjFileName,ProjFilePathName,EBSDFileName,EBSDFilePathName,CalFileName,CalFilePathName,filterVal,handles)
%Interpolates strain data to read into MTEX to enable combination with EBSD
%data. The output saves the DIC data into the MTEX structure. To set
%deform_switch you need to select whether the EBSD map you are reading in
%is a pre or post map. If ?reference? the original image positions are used
%and if ?deform? is selected then the deformed locations are used.
%
%This can only be used after using straincalc7 or above with deform_switch
%set to deform.

%
%Alex Forsey    Apr 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = waitbar(0,'straininterp');
if exist('ProjFileName','var')&&exist('ProjFilePathName','var')&&exist('EBSDFileName','var')&&exist('EBSDFilePathName','var')&&exist('CalFileName','var')&&exist('CalFilePathName','var')&&exist('filterVal','var')
    function_switch = 'function';
    load([EBSDFilePathName EBSDFileName],'grains','ebsd')
    load([ProjFilePathName ProjFileName])
    load([CalFilePathName CalFileName])
else
    function_switch = 'script';
    scale_switch = 'load';%'define';
    
    %outlier filter applied when calculating average strain in a grain from displacements.
    % remove data more than filterVal standard deviations from the mean
    filterVal = 2;%0;%9;%5;%0;%3;%0 = off
    %load strain tensor file and ebsd MTEX data
    [ProjFileName,ProjFilePathName] = uigetfile('*.mat','read DIC mat file');
    load([ProjFilePathName ProjFileName])
    [EBSDFileName,EBSDFilePathName] = uigetfile('*.mat','read ebsd mat file');
    load([EBSDFilePathName EBSDFileName],'grains','ebsd')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Set scale for DIC data
    switch scale_switch
        case 'load'
            [CalFileName,CalFilePathName] = uigetfile('*.mat','read calibration .mat file');
            load([CalFilePathName CalFileName])
            
        case 'define'
            %change this switch to use either the deformed position data or the
            %reference position data to compare to the ebsd
            deform_switch = 'deform';%'reference';%
            
            
            % DIC.mm2pixX = 9.25926e-005;%mm per pixel in x direction
            %directly from Da Vis this works only for HR5
            %LR5 = 9.25926e-005;
            DICmm2pixX = 3.7037e-005;%mm per pixel in x direction
            %DIC.mm2pixX = 9.25926e-005;%mm per pixel in x direction
            
            DICmm2pixY = DIC.mm2pixX;%mm per pixel in y direction
            %move in um backwards
            DICoriginX = 80;%0;%-60;%0;%0.13;%0.060;%x origin in pixels [corresponds to (0,0) in EBSD map]
            DICoriginY = 51.5;%0;%-40;%0;%.08;%0.030;%y origin in pixels [corresponds to (0,0) in EBSD map]
            
            save([EBSDFilePathName EBSDFileName(1:end-4) 'Calibration'],'DIC.mm2pixX','DIC.mm2pixY','DIC.originX','DIC.originY','deform_switch')
    end
    
end
DIC.mm2pixX = DICmm2pixX;
DIC.originX = DICoriginX;
DIC.mm2pixY = DICmm2pixY;
DIC.originY = DICoriginY;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EBSDposX = ebsd.x;
EBSDposY = ebsd.y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%apply scaling to DIC
switch deform_switch
    case 'reference'
        posX = reshape((DIC.Eposx*DIC.mm2pixX),[],1)*-1e3 + DIC.originX;
        posY = reshape((DIC.Eposy*DIC.mm2pixY),[],1)*-1e3 + DIC.originY;
        
    case 'deform'
        posX = reshape((DIC.Eposx_def*DIC.mm2pixX),[],1)*-1e3 + DIC.originX;
        posY = reshape((DIC.Eposy_def*DIC.mm2pixY),[],1)*-1e3 + DIC.originY;
end

%remove NaNs before interpolation
% Exx = reshape(squeeze(DIC.E(1,1,:)),[],1);
% Exy = reshape(squeeze(DIC.E(1,2,:)),[],1);
% Eyx = reshape(squeeze(DIC.E(2,1,:)),[],1);
% Eyy = reshape(squeeze(DIC.E(2,2,:)),[],1);
Exx = reshape(squeeze(DIC.E(1,1,:,:)),[],1);
Exy = reshape(squeeze(DIC.E(1,2,:,:)),[],1);
Eyx = reshape(squeeze(DIC.E(2,1,:,:)),[],1);
Eyy = reshape(squeeze(DIC.E(2,2,:,:)),[],1);



        posX(isnan(Eyy)) = [];posY(isnan(Eyy)) = [];
        Exx(isnan(Eyy)) = [];Exy(isnan(Eyy)) = [];Eyx(isnan(Eyy)) = [];Eyy(isnan(Eyy)) = [];
        Exx(isnan(posY)) = [];Exy(isnan(posY)) = [];Eyx(isnan(posY)) = [];Eyy(isnan(posY)) = [];
        posX(isnan(posY)) = [];posY(isnan(posY)) = [];

     

waitbar(1/10,H);
%interpolate tensor data
version -release;
Version = ans;
if str2double(Version(1:4))<=2011
    Fxx = TriScatteredInterp(posX,posY,Exx,'linear');waitbar(2/10,H);
    Fxy = TriScatteredInterp(posX,posY,Exy,'linear');waitbar(3/10,H);
    Fyx = TriScatteredInterp(posX,posY,Eyx,'linear');waitbar(4/10,H);
    Fyy = TriScatteredInterp(posX,posY,Eyy,'linear');waitbar(5/10,H);
else
    Fxx = scatteredInterpolant(posX,posY,Exx,'linear','nearest');waitbar(2/10,H);
    Fxy = scatteredInterpolant(posX,posY,Exy,'linear','nearest');waitbar(3/10,H);
    Fyx = scatteredInterpolant(posX,posY,Eyx,'linear','nearest');waitbar(4/10,H);
    Fyy = scatteredInterpolant(posX,posY,Eyy,'linear','nearest');waitbar(5/10,H);
end

% if str2double(Version(1:4))<=2011
%     Fxx = TriScatteredInterp(DIC.posX,DIC.posY,Exx,'linear');waitbar(2/10,H);
%     Fxy = TriScatteredInterp(DIC.posX,DIC.posY,Exy,'linear');waitbar(3/10,H);
%     Fyx = TriScatteredInterp(DIC.posX,DIC.posY,Eyx,'linear');waitbar(4/10,H);
%     Fyy = TriScatteredInterp(DIC.posX,DIC.posY,Eyy,'linear');waitbar(5/10,H);
% else
%     Fxx = scatteredInterpolant(DIC.posX,DIC.posY,Exx,'linear','nearest');waitbar(2/10,H);
%     Fxy = scatteredInterpolant(DIC.posX,DIC.posY,Exy,'linear','nearest');waitbar(3/10,H);
%     Fyx = scatteredInterpolant(DIC.posX,DIC.posY,Eyx,'linear','nearest');waitbar(4/10,H);
%     Fyy = scatteredInterpolant(DIC.posX,DIC.posY,Eyy,'linear','nearest');waitbar(5/10,H);
% end
Eint(1,1,:) = Fxx(EBSDposX,EBSDposY);waitbar(6/10,H);
Eint(1,2,:) = Fxy(EBSDposX,EBSDposY);waitbar(7/10,H);
Eint(2,1,:) = Fyx(EBSDposX,EBSDposY);waitbar(8/10,H);
Eint(2,2,:) = Fyy(EBSDposX,EBSDposY);waitbar(9/10,H);

%calculate principal strains and rotations
Ep = zeros(2,2,size(Eint,3),1);RotMat = Ep;Gmax = zeros(size(Eint,3),1);Theta = Gmax;
for jj = 1:size(Eint,3)
    [Ep(:,:,jj),Gmax(jj),Theta(jj),RotMat(:,:,jj),Nu(jj)] = principalStrains(squeeze(Eint(:,:,jj)));
end

setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','intoPlane');
setMTEXpref('yAxisDirection','south');

%save in MTEX format
ebsd.prop.Exx = squeeze(Eint(1,1,:));        %Interpolated Exx
ebsd.prop.Exy = squeeze(Eint(1,2,:));        %Interpolated Exy
ebsd.prop.Eyx = squeeze(Eint(2,1,:));        %Interpolated Eyx
ebsd.prop.Eyy = squeeze(Eint(2,2,:));        %Interpolated Eyy
ebsd.prop.EpMax = squeeze(Ep(1,1,:));        %Max principal strain
ebsd.prop.EpMin = squeeze(Ep(2,2,:));        %Min principal strain
ebsd.prop.Rotation = squeeze(RotMat(1,2,:)); %Rigid body rotation
ebsd.prop.Gmax = Gmax;              %Principal shear strain
ebsd.prop.Theta = Theta;            %Principal strain angle
ebsd.prop.Poisson = Nu;        %Poisson's ratio


%crop data to extremes of the DIC strain data
% % xcrop=[max(DIC.posX(:)),min(DIC.posX(:))];
% % ycrop=[max(DIC.posY(:)),min(DIC.posY(:))];
% % ebsd=ebsd(ebsd.x>xcrop(1) & ebsd.x<xcrop(2) & ebsd.y>ycrop(1) & ebsd.y<ycrop(2));
% % [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',5*degree);
% % [grains,DIC] = meanGrainStrain(ebsd,grains,DIC,deform_switch,filterVal);
waitbar(1,H);
close(H)
cropLoop = 'no';
while strcmp(cropLoop,'no')
    %figure to plot from
    figure
    h = gca;
    set(h,'position',[377 17 710 533])
    
    plot(ebsd,ebsd.Exx)
    caxis([0 0.3])
    hold on
    plot(grains.boundary,'linewidth',2)
    mtexColorMap hot
    hold off
    title('Click opposing corners to crop plots')
    
    %crop data for plotting
    [x,y]= ginput(1);
    hold on
    plot(x,y,'og','LineWidth',2)
    [x(2),y(2)]= ginput(1);
    plot(x(2),y(2),'og','LineWidth',2)
    hold off
    close
    x=sort(x);y=sort(y);
%     ebsd=ebsd(ebsd.x>x(1) & ebsd.x<x(2) & ebsd.y>y(1) & ebsd.y<y(2));
    ebsd2=ebsd(ebsd.x>x(1) & ebsd.x<x(2) & ebsd.y>y(1) & ebsd.y<y(2));
    [grains2,ebsd2.grainId,ebsd2.mis2mean] = calcGrains(ebsd2,'angle',5*degree);
    figure
    h = gca;
    set(h,'position',[377 17 710 533])
    plot(ebsd2,ebsd2.Exx)
    
%     plot(ebsd,ebsd.Exx)
    hold on
    plot(grains2.boundary,'linewidth',2)
    caxis([0 0.3])
    mtexColorMap hot
    hold off
    title('Exx')
    
    cropLoop = questdlg('Accept crop?', 'Crop control', 'yes', 'no', 'yes');
    close
end

if exist('grains2','var')
    grains = grains2;
    ebsd = ebsd2;
    clear grains2 ebsd2
end

if isfield(DIC,'raw')
    
    [grains,DIC] = meanGrainStrain(ebsd,grains,DIC,deform_switch,filterVal);
else
    grains = grainStrain(ebsd,grains,H);
end
            
save([ProjFilePathName ProjFileName(1:end-4) '_ebsdDIC_std' num2str(filterVal)],'ebsd','grains','DIC')


end

function [e1,Gmax,Theta,W,Nu] = principalStrains(e)
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

%calculate Poisson's ratio
% Nu = min(abs(ee([1,4]))/max(abs(ee([1,4]))));
Nu = -ee(2,2)/ee(1,1);

end

function [grains,DIC] = meanGrainStrain(ebsd,grains,DIC,deform_switch,filterVal)
%interpolate grain reference number for original displacement vectors using
%nearest approach

%remove nans
DIC.raw.x = DIC.raw.x(~isnan(DIC.raw.u));
DIC.raw.y = DIC.raw.y(~isnan(DIC.raw.u));
DIC.raw.v = DIC.raw.v(~isnan(DIC.raw.u));
DIC.raw.u = DIC.raw.u(~isnan(DIC.raw.u));

%apply scaing
DIC.raw.x = DIC.raw.x*DIC.mm2pixX*-1e3 + DIC.originX;
DIC.raw.y = DIC.raw.y*DIC.mm2pixY*-1e3 + DIC.originY;
DIC.raw.u = DIC.raw.u*DIC.mm2pixX*-1e3;
DIC.raw.v = DIC.raw.v*DIC.mm2pixY*-1e3;

%interpolate grain Id
version -release;
Version = ans;
if str2double(Version(1:4))<=2011
    FgrainId = TriScatteredInterp(ebsd.x,ebsd.y,ebsd.grainId,'nearest');
else
    FgrainId = scatteredInterpolant(ebsd.x,ebsd.y,ebsd.grainId,'nearest','nearest');
end
switch deform_switch
    case 'reference'
        DIC.raw.grainId = FgrainId(DIC.raw.x,DIC.raw.y);
    case 'deform'
        DIC.raw.grainId = FgrainId(DIC.raw.x + DIC.raw.u,DIC.raw.y + DIC.raw.v);
end

grains.prop.Exx=zeros(size(grains));
grains.prop.Eyy=zeros(size(grains));
grains.prop.Exy=zeros(size(grains));
grains.prop.Eyx=zeros(size(grains));
grains.prop.EpMax=zeros(size(grains));
grains.prop.EpMin=zeros(size(grains));
grains.prop.Rotation=zeros(size(grains));
grains.prop.Gmax=zeros(size(grains));
grains.prop.Theta=zeros(size(grains));
grains.prop.Poisson=zeros(size(grains));

%for each grain calculate the average strain from the displacement vectors
%for that grain
H = waitbar(0,'straininterp');
for i = 1:length(grains)
    waitbar(i/length(grains),H)
    %identify vectors whithin active grain
    singleGrainCrop = DIC.raw.grainId==i;
    %fit a plane to the displacements to calculate strain tensor
    if max(singleGrainCrop)~=0
        [E,Coeffs] = etensor3(DIC.raw.x(singleGrainCrop),DIC.raw.y(singleGrainCrop),DIC.raw.u(singleGrainCrop),DIC.raw.v(singleGrainCrop));
        %if required remove outliers and fit again (filterVal of 0 turns it off)
        if filterVal ~= 0 && max(~isnan(E(:)))
            %create ideal vector field for measured strain tensor
            [xfield,yfield] = flowmaker3(Coeffs,DIC.raw.x,DIC.raw.y);
            %calculate difference between measured and ideal data
            xdiff = DIC.raw.u-xfield;
            ydiff = DIC.raw.v-yfield;
            %create outlier filter
            stdxdiff = std(xdiff);stdydiff = std(ydiff);
            
            outlierCropx1 = xdiff>+filterVal*stdxdiff;
            outlierCropx2 = xdiff<-filterVal*stdxdiff;
            outlierCropy1 = ydiff>+filterVal*stdydiff;
            outlierCropy2 = ydiff<-filterVal*stdydiff;
            %combine crops
            combCrop = singleGrainCrop&~(outlierCropx1|outlierCropx2|outlierCropy1|outlierCropy2);
            %recalculate strain tensor
            [E,~] = etensor3(DIC.raw.x(combCrop),DIC.raw.y(combCrop),DIC.raw.u(combCrop),DIC.raw.v(combCrop));
        end
    else
        E = nan(2);
    end
    
    %calculate principal strains
    [Ep,Gmax,Theta,W,Nu] = principalStrains(E);
    
    %save into MTex structure
    
    grains(i).Exx = E(1,1);
    grains(i).Eyy = E(2,2);
    grains(i).Exy = E(1,2);
    grains(i).Eyx = E(2,1);
    grains(i).EpMax = Ep(1,1);
    grains(i).EpMin = Ep(2,2);
    grains(i).Rotation = W(1,2);
    grains(i).Gmax = Gmax;
    grains(i).Theta = Theta;
    grains(i).Poisson = Nu;
end
close(H)
end

function [E,Coeffs] = etensor3(X,Y,U,V)
X = X(:); % Make X a column vector
Y = Y(:); % Make Y a column vector
U = U(:); % Make Z a column vector
V = V(:); % Make Z a column vector

origLengthX = length(X);

X(isnan(U)) = [];
Y(isnan(U)) = [];
V(isnan(U)) = [];
U(isnan(U)) = [];

grainFilter = 'no';
switch grainFilter
    case 'yes'
        %don't calculate if there are insufficient data points in either direction
        %or if over half of the window is unpopulated
        if length(unique(X))>=2 && length(unique(Y))>=2 && length(X)/origLengthX > 0.5
            warning off MATLAB:nearlySingularMatrix;
            Const = ones(size(X)); % Vector of ones for constant term
            CoefficientsU = [X Y Const]\U; % Find the coefficients for x displacement
            CoefficientsV = [X Y Const]\V; % Find the coefficients for y displacement
            
            %remove low-confidence points
            [~, msgidlast] = lastwarn;
            if strcmp(msgidlast,'MATLAB:nearlySingularMatrix')
                lastwarn('')
                E = nan(2);
                Coeffs = {NaN,NaN};
            else
                E=[CoefficientsU(1),CoefficientsU(2);CoefficientsV(1),CoefficientsV(2)];
                Coeffs = {CoefficientsU,CoefficientsV};
            end
            warning on MATLAB:nearlySingularMatrix;
            
        else
            E = nan(2);
            Coeffs = {NaN,NaN};
        end
    case 'no' % save all grain average data regardless of quality
        warning off MATLAB:nearlySingularMatrix;
        Const = ones(size(X)); % Vector of ones for constant term
        CoefficientsU = [X Y Const]\U; % Find the coefficients for x displacement
        CoefficientsV = [X Y Const]\V; % Find the coefficients for y displacement
        E=[CoefficientsU(1),CoefficientsU(2);CoefficientsV(1),CoefficientsV(2)];
        Coeffs = {CoefficientsU,CoefficientsV};
        warning on MATLAB:nearlySingularMatrix;
end
end

function [xfield,yfield] = flowmaker3(Coeffs,X,Y)
%Interpolate fitted data parameter to given points
xfield = Coeffs{1}(1) .* X + Coeffs{1}(2) .* Y + Coeffs{1}(3);
yfield = Coeffs{2}(1) .* X + Coeffs{2}(2) .* Y + Coeffs{2}(3);
end

function grains = grainStrain(ebsd,grains,H)
%calculate average strain in each grain from MTEX
gId=ebsd.grainId;
grains.prop.Exx=zeros(size(grains));
grains.prop.Eyy=zeros(size(grains));
grains.prop.Exy=zeros(size(grains));
grains.prop.Eyx=zeros(size(grains));
grains.prop.EpMax=zeros(size(grains));
grains.prop.EpMin=zeros(size(grains));
grains.prop.Rotation=zeros(size(grains));
grains.prop.Gmax=zeros(size(grains));
grains.prop.Theta=zeros(size(grains));
grains.prop.Poisson=zeros(size(grains));

for n=1:length(grains)
    waitbar(n/length(grains),H)
    ebsdGpos=find(gId==n);
    
    grains(n).Exx=mean( ebsd( ebsdGpos ).Exx );
    grains(n).Eyy=mean( ebsd( ebsdGpos ).Eyy );
    grains(n).Exy=mean( ebsd( ebsdGpos ).Exy );
    grains(n).Eyx=mean( ebsd( ebsdGpos ).Eyx );
    grains(n).EpMax=mean( ebsd( ebsdGpos ).EpMax );
    grains(n).EpMin=mean( ebsd( ebsdGpos ).EpMin );
    grains(n).Rotation=mean( ebsd( ebsdGpos ).Rotation );
    grains(n).Gmax=mean( ebsd( ebsdGpos ).Gmax );
    grains(n).Theta=mean( ebsd( ebsdGpos ).Theta );
    grains(n).Poisson=mean( ebsd( ebsdGpos ).Poisson );
end
end

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
        cropLoop = 'yes';
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