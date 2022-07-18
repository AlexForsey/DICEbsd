function DIC = straincalc6_2GUI(DIC,handles)
%Calculates strain using a moving strain window on a single exported (.txt)
%DIC data. Four parts of the strain tensor are saved and the true position
%of strain is also included.
%
%Output is four dimensional matrix E.
%Dimension 1 is x position
%Dimension 2 is y position
%Dimensions 3 and 4 describe the strain tensor.
% ie E(e11,e12;
%      e21,e22)
%
%
%
%Alex Forsey April 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Set strain window size
% SW = 2;%increase this value to increase smoothing (default = 1 and gives strain calculated from 2x2 vectors)
% DIC.SW = [SW,SW];%change here to smooth different ammounts in different directions by definins SW in y [SW(1)] and x [SW(2)]

%save deformed data point locations (for MTEX)?
deform_switch = 'deform';%'reference only';%
principal_switch = 'off';%'calculate';%
save_displacements = 'on';%'off';%

%Select file
% [DIC.ProjFileName,DIC.ProjFilePathName] = uigetfile('*.txt','read text file','MultiSelect', 'on');

%calculate region used for strain calculation
% window =  input('Input window size (pixels)');
% step = input('Input step size (pixels)');
% DIC.SWpix = window + SW*step;%the length of one side of the strain window


if iscell(DIC.ProjFileName)
    totalFiles = length(DIC.ProjFileName);
else
    totalFiles = 1;
end
overlapSwitch = get(handles.overlapSW,'Value');

for  ijk = 1:totalFiles
    %Read in txt file
    if totalFiles == 1
        fname = DIC.ProjFileName;
    else
        fname = cell2mat(DIC.ProjFileName(ijk));
    end
    [x,y,u,v] = DICtxtReadSub([DIC.ProjFilePathName,fname]);
    v(u==0) = nan;
    u(u==0) = nan;
    
    %Strain calc
%     DIC.Eposy=zeros(size(x,1)-DIC.SW(1),size(x,2)-DIC.SW(2));
%     DIC.Eposx=zeros(size(x,1)-DIC.SW(1),size(x,2)-DIC.SW(2));
    H = waitbar(0,['Straincalc    ...' fname]);wb = 1;
    switch overlapSwitch
        case 1
            iLim = size(x,1)-DIC.SW(1);
            jLim = size(x,2)-DIC.SW(2);
        case 0
            iLim = floor((size(x,1)-DIC.SW(1))/DIC.SW(1));
            jLim = floor((size(x,2)-DIC.SW(2))/DIC.SW(2));
%             iLim = floor((size(x,1))/DIC.SW(1));
%             jLim = floor((size(x,2))/DIC.SW(2));
    end
    
    for i = 1:iLim
        for j = 1:jLim
            waitbar(wb/(iLim*jLim),H);wb = wb+1;
            switch overlapSwitch
                case 0
                    %crop to strain window
                    xCrop = (i*DIC.SW(1)-DIC.SW(1)+1:i*DIC.SW(1)+1);
                    yCrop = (j*DIC.SW(2)-DIC.SW(2)+1:j*DIC.SW(2)+1);
                case 1
                    %crop to strain window
                    xCrop = (i:i+DIC.SW(1));
                    yCrop = (j:j+DIC.SW(2));
                    
            end
            
            %calculate strain tensor
            [DIC.E(:,:,i,j),Coeffs] = etensor3(x(xCrop,yCrop),y(xCrop,yCrop),u(xCrop,yCrop),v(xCrop,yCrop));
            
            %Calculate principal strain data
            switch principal_switch
                case 'calculate'
                    [DIC.Ep(:,:,i,j),DIC.Gmax(i,j),DIC.Theta(i,j),DIC.RotMat(:,:,i,j),DIC.Poisson(i,j)] = principalStrains(squeeze(DIC.E(:,:,i,j)));
            end
            
            %calculate position of strain value
            DIC.Eposx(i,j) = mean(mean(x(xCrop,yCrop)));
            DIC.Eposy(i,j) = mean(mean(y(xCrop,yCrop)));
        end
    end
    
    %plot strain
    axes(handles.axes1)
%     handles.axes1 = imagesc(DIC.Eposx(1,:),DIC.Eposy(:,1),squeeze(DIC.E(1,1,:,:))*100);
imagesc(DIC.Eposx(1,:),DIC.Eposy(:,1),squeeze(DIC.E(1,1,:,:))*100);
    axis image
    colorbar
    Cax(1) = str2double(get(handles.caxmin,'string'));
    Cax(2) = str2double(get(handles.caxmax,'string'));
    if ~isnan(mean(Cax(:)))&&Cax(1)<Cax(2) 
        caxis(Cax)
    else 
        Cax = caxis;
        
        set(handles.caxmin,'string',num2str(Cax(1)))
        set(handles.caxmax,'string',num2str(Cax(2)))
    end
    
    switch deform_switch
        case 'deform'
%             switch overlapSwitch
%                 case 1
%                     DIC.Eposx_def = DIC.Eposx+u(ceil(DIC.SW(2)/2):end-ceil(DIC.SW(2)/2),ceil(DIC.SW(2)/2):end-ceil(DIC.SW(2)/2));
%                     DIC.Eposy_def = DIC.Eposy+v(ceil(DIC.SW(1)/2):end-ceil(DIC.SW(1)/2),ceil(DIC.SW(1)/2):end-ceil(DIC.SW(1)/2));
%     
%                 case 0
                    Fu = TriScatteredInterp(x(:),y(:),u(:));
                    Fv = TriScatteredInterp(x(:),y(:),v(:));
                    
                    DIC.Eposx_def = DIC.Eposx+Fu(DIC.Eposx,DIC.Eposy);
                    DIC.Eposy_def = DIC.Eposy+Fv(DIC.Eposx,DIC.Eposy);
    
%             end
    end
    switch save_displacements
        case 'on'
            DIC.raw.x = x(:);
            DIC.raw.y = y(:);
            DIC.raw.u = u(:);
            DIC.raw.v = v(:);
    end
    DIC.fname = ['Strain_SW_' num2str(DIC.SW(1)) '_' num2str(DIC.SW(2)) '_' fname(1:end-4)];
    switch overlapSwitch
        case 1
            save([DIC.ProjFilePathName, DIC.fname],'DIC')
        case 0
            save([DIC.ProjFilePathName, DIC.fname,'_sparse'],'DIC')
     
    end
end
close(H)
end


function [x,y,u,v] = DICtxtReadSub(fname)
%Reads in .txt files and orders them into matrices for processing
% data = dlmread(fname,'\t',1,0);;
data = importdata(fname);
data = data.data;

%sort out orientation of data
[nC,nIA,nIC] = (unique(data(:,1)));
[mC,mIA,mIC] = (unique(data(:,2)));
n = length(nC);m = length(mC);
x = nan(m,n);y = x;u = x;v = x;

for o = 1:length(data)
    x(mIC(o),nIC(o)) = data(o,1);
    y(mIC(o),nIC(o)) = data(o,2);
    u(mIC(o),nIC(o)) = data(o,3);
    v(mIC(o),nIC(o)) = data(o,4);
end



% x = transpose(reshape(data(:,1),m,n));
% y = transpose(reshape(data(:,2),m,n));
% u = transpose(reshape(data(:,3),m,n));
% v = transpose(reshape(data(:,4),m,n));
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

end

function [e1,Gmax,Theta,W,Poisson] = principalStrains(e)
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
Poisson = min(abs(ee([1,4]))/max(abs(ee([1,4]))));
end

function [xfield,yfield] = flowmaker3(Coeffs,X,Y)
%Interpolate fitted data parameter to given points
xfield = Coeffs{1}(1) .* X + Coeffs{1}(2) .* Y + Coeffs{1}(3);
yfield = Coeffs{2}(1) .* X + Coeffs{2}(2) .* Y + Coeffs{2}(3);
end