% function DICEscale(MTEX)
function varargout = DICEscale(varargin)
% DICESCALE MATLAB code for DICEscale.fig
%      DICESCALE, by itself, creates a new DICESCALE or raises the existing
%      singleton*.
%
%      H = DICESCALE returns the handle to a new DICESCALE or the handle to
%      the existing singleton*.
%
%      DICESCALE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DICESCALE.M with the given input arguments.
%
%      DICESCALE('Property','Value',...) creates a new DICESCALE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DICEscale_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DICEscale_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DICEscale

% Last Modified by GUIDE v2.5 23-May-2022 09:13:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DICEscale_OpeningFcn, ...
                   'gui_OutputFcn',  @DICEscale_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%save input data to gui 


% --- Executes just before DICEscale is made visible.
function DICEscale_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DICEscale (see VARARGIN)

% Choose default command line output for DICEscale
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DICEscale wait for user response (see UIRESUME)
% uiwait(handles.DiceScale);
hDICE = findobj('Tag','dice');
MTEX = getappdata(hDICE,'MTEX');
load([MTEX.EBSD.ProjFilePathName MTEX.EBSD.ProjFileName],'grains','ebsd')
load([MTEX.DIC.ProjFilePathName MTEX.DIC.ProjFileName],'DIC')

setappdata(handles.output,'MTEX',MTEX)
setappdata(handles.output,'DIC',DIC)
setappdata(handles.output,'ebsd',ebsd)
setappdata(handles.output,'grains',grains)


% --- Outputs from this function are returned to the command line.
function varargout = DICEscale_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function scaleFactor_Callback(hObject, eventdata, handles)
% hObject    handle to scaleFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scaleFactor as text
%        str2double(get(hObject,'String')) returns contents of scaleFactor as a double


% --- Executes during object creation, after setting all properties.
function scaleFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scaleFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xOrigin_Callback(hObject, eventdata, handles)
% hObject    handle to xOrigin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xOrigin as text
%        str2double(get(hObject,'String')) returns contents of xOrigin as a double


% --- Executes during object creation, after setting all properties.
function xOrigin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xOrigin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yOrigin_Callback(hObject, eventdata, handles)
% hObject    handle to yOrigin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yOrigin as text
%        str2double(get(hObject,'String')) returns contents of yOrigin as a double


% --- Executes during object creation, after setting all properties.
function yOrigin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yOrigin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in deformData.
function deformData_Callback(hObject, eventdata, handles)
% hObject    handle to deformData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of deformData



function calFname_Callback(hObject, eventdata, handles)
% hObject    handle to calFname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calFname as text
%        str2double(get(hObject,'String')) returns contents of calFname as a double


% --- Executes during object creation, after setting all properties.
function calFname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calFname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in testButton.
function testButton_Callback(hObject, eventdata, handles)
% hObject    handle to testButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get plot data
DIC = getappdata(handles.output,'DIC');
ebsd = getappdata(handles.output,'ebsd');
grains = getappdata(handles.output,'grains');
DICmm2pixX = str2double(get(handles.scaleFactor,'String'))*1e3;
DICmm2pixY = DICmm2pixX;
DICoriginX = str2double(get(handles.xOrigin,'String'));
DICoriginY = str2double(get(handles.yOrigin,'String'));
calFname = get(handles.calFname,'String');
defSw = get(handles.deformData,'Value');
plotPropStr = get(handles.plotProp,'String');
plotPropVal = get(handles.plotProp,'Value');
plotProp = cell2mat(plotPropStr(plotPropVal));

if defSw == 1
    deform_switch = 'deform';
else
    deform_switch = 'reference';
end

% affine
a2 = str2double(get(handles.a2,'String'));
a3 = str2double(get(handles.a3,'String'));
a4 = str2double(get(handles.a4,'String'));
a5 = str2double(get(handles.a5,'String'));
% a6 = 0;
% a7 = 0;
% a8 = 1;

TransMat = [a2,a3;a4,a5];
% TransMat3 = [a2,a3,a6;a4,a5,a7;a6,a7,a8];

EBSDcoords = TransMat*[ebsd.prop.x';ebsd.prop.y'];
% EBSDcoords = TransMat3*[ebsd.prop.x';ebsd.prop.y';ones(size(ebsd.prop.x'))];
% Grainscoords = TransMat*[grains.x';grains.y'];

% figure
% plot(coords(:,1),coords(:,2),'.')

ebsd.prop.x = EBSDcoords(1,:)';
ebsd.prop.y = EBSDcoords(2,:)';

grains = calcGrains(ebsd,'angle',15*degree);
    
% plotting
plotHandle = findobj('Tag','scalesAxes');
axes(plotHandle)

% Flip map switch
if ~get(handles.rvseX,'Value')
    DICmm2pixX = DICmm2pixX*(-1);
end
if ~get(handles.rvseY,'Value')
    DICmm2pixY = DICmm2pixY*(-1);
end

figure(1)
hold off
plot(grains.boundary)
hold on

%select strain value to plot
if strcmp(plotProp,'Exx')
    Eplot = squeeze(DIC.E(1,1,:,:));
elseif strcmp(plotProp,'Exy')
    Eplot = squeeze(DIC.E(2,1,:,:));
elseif strcmp(plotProp,'Eyx')
    Eplot = squeeze(DIC.E(1,2,:,:));
% elseif strcmp(plotProp,'Eyy')
else 
    Eplot = squeeze(DIC.E(2,2,:,:));
end

% if strcmp(deform_switch,'deform')
%     %interpolate deformed data to grid to plot using imagesc
%         Fdef = scatteredInterpolant(DIC.Eposx_def(~isnan((DIC.Eposx_def))),DIC.Eposy_def(~isnan((DIC.Eposx_def))),Eplot(~isnan((DIC.Eposx_def))));
%         Edef = Fdef(DIC.Eposx,DIC.Eposy);
%     imagesc(((DIC.Eposx(1,:).*DICmm2pixX)+DICoriginX),((DIC.Eposy(:,1).*DICmm2pixY)+DICoriginY),Edef)
% else
% %     imagesc((DIC.Eposx(1,:).*DICmm2pixX)+DICoriginX,(DIC.Eposy(:,1).*DICmm2pixY)+DICoriginY,Eplot)
% end

imagesc((DIC.image.x.*DICmm2pixX)+DICoriginX,(DIC.image.y.*DICmm2pixY)+DICoriginY,DIC.image.I)


% caxis([str2double(get(handles.LcAxis,'String')),str2double(get(handles.UcAxis,'String'))])
plot(grains.boundary)
hold off

% Flip axes and units to agree with MTex later on - remove at some point for consistency 
DICmm2pixX = DICmm2pixX*(-1e-3);
DICmm2pixY = DICmm2pixY*(-1e-3);

save(calFname,'DICoriginX','DICoriginY','DICmm2pixX','DICmm2pixY','deform_switch')

% if TransMat ~= eye(2)
    save('affineEBSD.mat','ebsd','grains')
% end


% if 
% save(calFname,'DICoriginX','DICoriginY','DICmm2pixX','DICmm2pixY','deform_switch')



function UcAxis_Callback(hObject, eventdata, handles)
% hObject    handle to UcAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UcAxis as text
%        str2double(get(hObject,'String')) returns contents of UcAxis as a double


% --- Executes during object creation, after setting all properties.
function UcAxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UcAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LcAxis_Callback(hObject, eventdata, handles)
% hObject    handle to LcAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LcAxis as text
%        str2double(get(hObject,'String')) returns contents of LcAxis as a double


% --- Executes during object creation, after setting all properties.
function LcAxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LcAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rvseY.
function rvseY_Callback(hObject, eventdata, handles)
% hObject    handle to rvseY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rvseY


% --- Executes on button press in rvseX.
function rvseX_Callback(hObject, eventdata, handles)
% hObject    handle to rvseX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rvseX


% --- Executes on selection change in plotProp.
function plotProp_Callback(hObject, eventdata, handles)
% hObject    handle to plotProp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plotProp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotProp


% --- Executes during object creation, after setting all properties.
function plotProp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotProp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a2_Callback(hObject, eventdata, handles)
% hObject    handle to a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a2 as text
%        str2double(get(hObject,'String')) returns contents of a2 as a double


% --- Executes during object creation, after setting all properties.
function a2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a3_Callback(hObject, eventdata, handles)
% hObject    handle to a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a3 as text
%        str2double(get(hObject,'String')) returns contents of a3 as a double


% --- Executes during object creation, after setting all properties.
function a3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a4_Callback(hObject, eventdata, handles)
% hObject    handle to a4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a4 as text
%        str2double(get(hObject,'String')) returns contents of a4 as a double


% --- Executes during object creation, after setting all properties.
function a4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a5_Callback(hObject, eventdata, handles)
% hObject    handle to a5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a5 as text
%        str2double(get(hObject,'String')) returns contents of a5 as a double


% --- Executes during object creation, after setting all properties.
function a5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
