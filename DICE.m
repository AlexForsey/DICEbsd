function varargout = DICE(varargin)
% DICE MATLAB code for DICE.fig
%      DICE, by itself, creates a new DICE or raises the existing
%      singleton*.
%
%      H = DICE returns the handle to a new DICE or the handle to
%      the existing singleton*.
%
%      DICE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DICE.M with the given input arguments.
%
%      DICE('Property','Value',...) creates a new DICE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DICE_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DICE_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DICE

% Last Modified by GUIDE v2.5 02-Mar-2018 12:23:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DICE_OpeningFcn, ...
                   'gui_OutputFcn',  @DICE_OutputFcn, ...
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


% --- Executes just before DICE is made visible.
function DICE_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DICE (see VARARGIN)

% Choose default command line output for DICE
handles.output = hObject;
setappdata(handles.output,'DIC',[])
setappdata(handles.output,'MTEX',[])
setappdata(handles.output,'plotting',[])

% Update handles structure
guidata(hObject, handles);
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','inToPlane');
setMTEXpref('showMicronBar','off')
setMTEXpref('showCoordinates','on')

% UIWAIT makes DICE wait for user response (see UIRESUME)
% uiwait(handles.dice);


% --- Outputs from this function are returned to the command line.
function varargout = DICE_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function resultsFile_Callback(hObject, eventdata, handles)
% hObject    handle to resultsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resultsFile as text
%        str2double(get(hObject,'String')) returns contents of resultsFile as a double


% --- Executes during object creation, after setting all properties.
function resultsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resultsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browseResults.
function browseResults_Callback(hObject, eventdata, handles)
% hObject    handle to browseResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[plotting.ProjFileName,plotting.ProjFilePathName] = uigetfile('*.mat','read text file');
setappdata(handles.output,'plotting',plotting)
set(handles.resultsFile,'String',plotting.ProjFileName)


function caxmin_Callback(hObject, eventdata, handles)
% hObject    handle to caxmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of caxmin as text
%        str2double(get(hObject,'String')) returns contents of caxmin as a double
Cax(1) = str2double(get(handles.caxmin,'string'));
Cax(2) = str2double(get(handles.caxmax,'string'));
caxis(gca,Cax)

% --- Executes during object creation, after setting all properties.
function caxmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to caxmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function caxmax_Callback(hObject, eventdata, handles)
% hObject    handle to caxmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of caxmax as text
%        str2double(get(hObject,'String')) returns contents of caxmax as a double
Cax(1) = str2double(get(handles.caxmin,'string'));
Cax(2) = str2double(get(handles.caxmax,'string'));
caxis(gca,Cax)

% --- Executes during object creation, after setting all properties.
function caxmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to caxmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in autoScale.
function autoScale_Callback(hObject, eventdata, handles)
% hObject    handle to autoScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
caxis(gca,'auto')
Cax = caxis;

set(handles.caxmin,'string',num2str(Cax(1)))
set(handles.caxmax,'string',num2str(Cax(2)))


function p2_t_Callback(hObject, eventdata, handles)
% hObject    handle to p2_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p2_t as text
%        str2double(get(hObject,'String')) returns contents of p2_t as a double


% --- Executes during object creation, after setting all properties.
function p2_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p2_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p3_t_Callback(hObject, eventdata, handles)
% hObject    handle to p3_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p3_t as text
%        str2double(get(hObject,'String')) returns contents of p3_t as a double


% --- Executes during object creation, after setting all properties.
function p3_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p3_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p4_t_Callback(hObject, eventdata, handles)
% hObject    handle to p4_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p4_t as text
%        str2double(get(hObject,'String')) returns contents of p4_t as a double


% --- Executes during object creation, after setting all properties.
function p4_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p4_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p1_t_Callback(hObject, eventdata, handles)
% hObject    handle to p1_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p1_t as text
%        str2double(get(hObject,'String')) returns contents of p1_t as a double


% --- Executes during object creation, after setting all properties.
function p1_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p1_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plotType.
function plotType_Callback(hObject, eventdata, handles)
% hObject    handle to plotType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plotType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotType


% --- Executes during object creation, after setting all properties.
function plotType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotButton.
function plotButton_Callback(hObject, eventdata, handles)
% hObject    handle to plotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%DICEplot(handles,hObject);
%DICEplot(handles);%,hObject);
%DICEplot2(handles);%,hObject);
DICEplot3(handles);%,hObject);


cmax = str2double(get(handles.caxmax,'string'));
cmin = str2double(get(handles.caxmin,'string'));
if isempty(cmax)||isempty(cmin)||isnan(cmax)||isnan(cmin)
    caxis auto
    CAxis = caxis;
    set(handles.caxmax,'string',max(CAxis));
    set(handles.caxmin,'string',min(CAxis));
else
    caxis([cmin cmax])
end

% --- Executes on selection change in colorMap_m_pop.
function colorMap_m_pop_Callback(hObject, eventdata, handles)
% hObject    handle to colorMap_m_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns colorMap_m_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from colorMap_m_pop
colcs=get(handles.colorMap_m_pop,'string');choi=get(handles.colorMap_m_pop,'value');
%colcs=get(handles.colorMap_m,'string');choi=get(handles.colorMap_m,'value');
colc=colcs{choi};if choi==6; colc='LaboTeX';end

mtexColorMap(colc)

% --- Executes during object creation, after setting all properties.
function colorMap_m_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorMap_m_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plotMod.
function plotMod_Callback(hObject, eventdata, handles)
% hObject    handle to plotMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plotMod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotMod
varString = get(handles.plotMod,'string');
varVal = get(handles.plotMod,'Value');
varSwitch = cell2mat(varString(varVal));
if strcmp(varSwitch,'Scatter')
    set(handles.xAxisVar,'visible','on')
    set(handles.colorMap_m_pop,'visible','off')
    set(handles.text28,'string','y axis')
else
    set(handles.xAxisVar,'visible','off') 
    set(handles.colorMap_m_pop,'visible','on')
    set(handles.text28,'string','Colour map')
end

% --- Executes during object creation, after setting all properties.
function plotMod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in docPlot.
function docPlot_Callback(hObject, eventdata, handles)
% hObject    handle to docPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of docPlot
plotPosSwitch = get(handles.docPlot,'Value');
switch plotPosSwitch
    case 1
    set(handles.p1_t,'visible','off')
    set(handles.p2_t,'visible','off')
    set(handles.p3_t,'visible','off')
    set(handles.p4_t,'visible','off')
    set(handles.text23,'visible','off')
    case 0
    set(handles.p1_t,'visible','on')
    set(handles.p2_t,'visible','on')
    set(handles.p3_t,'visible','on')
    set(handles.p4_t,'visible','on')
    set(handles.text23,'visible','on')
end

% --- Executes on selection change in xAxisVar.
function xAxisVar_Callback(hObject, eventdata, handles)
% hObject    handle to xAxisVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xAxisVar contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xAxisVar


% --- Executes during object creation, after setting all properties.
function xAxisVar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xAxisVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MTEXname_Callback(hObject, eventdata, handles)
% hObject    handle to MTEXname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MTEXname as text
%        str2double(get(hObject,'String')) returns contents of MTEXname as a double


% --- Executes during object creation, after setting all properties.
function MTEXname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MTEXname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in createCALmat.
function createCALmat_Callback(hObject, eventdata, handles)
% hObject    handle to createCALmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MTEX = getappdata(handles.output,'MTEX');
% load([MTEX.EBSD.ProjFilePathName MTEX.EBSD.ProjFileName],'grains','ebsd')
% load([MTEX.DIC.ProjFilePathName MTEX.DIC.ProjFileName],'DIC')
% DICEscale(MTEX)DIC,ebsd,grains)
DICEscale(handles)

% DICEscale(handles)

% --- Executes on button press in browseCALmat.
function browseCALmat_Callback(hObject, eventdata, handles)
% hObject    handle to browseCALmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MTEX = getappdata(handles.output,'MTEX');
[MTEX.Cal.ProjFileName,MTEX.Cal.ProjFilePathName] = uigetfile('*.mat','read text file');
setappdata(handles.output,'MTEX',MTEX)
set(handles.calName,'String',MTEX.Cal.ProjFileName)


function calName_Callback(hObject, eventdata, handles)
% hObject    handle to calName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calName as text
%        str2double(get(hObject,'String')) returns contents of calName as a double


% --- Executes during object creation, after setting all properties.
function calName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browseMTEXmat.
function browseMTEXmat_Callback(hObject, eventdata, handles)
% hObject    handle to browseMTEXmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MTEX = getappdata(handles.output,'MTEX');
[MTEX.EBSD.ProjFileName,MTEX.EBSD.ProjFilePathName] = uigetfile('*.mat','read text file');
setappdata(handles.output,'MTEX',MTEX)
set(handles.MTEXname,'String',MTEX.EBSD.ProjFileName)

% --- Executes on button press in browseDICmat.
function browseDICmat_Callback(hObject, eventdata, handles)
% hObject    handle to browseDICmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MTEX = getappdata(handles.output,'MTEX');
[MTEX.DIC.ProjFileName,MTEX.DIC.ProjFilePathName] = uigetfile('*.mat','read text file');
setappdata(handles.output,'MTEX',MTEX)
set(handles.loadDICmat,'String',MTEX.DIC.ProjFileName)
 
load([MTEX.DIC.ProjFilePathName MTEX.DIC.ProjFileName])

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

function loadDICmat_Callback(hObject, eventdata, handles)
% hObject    handle to loadDICmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadDICmat as text
%        str2double(get(hObject,'String')) returns contents of loadDICmat as a double


% --- Executes during object creation, after setting all properties.
function loadDICmat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadDICmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filterVal_Callback(hObject, eventdata, handles)
% hObject    handle to filterVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filterVal as text
%        str2double(get(hObject,'String')) returns contents of filterVal as a double


% --- Executes during object creation, after setting all properties.
function filterVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in runInterp.
function runInterp_Callback(hObject, eventdata, handles)
% hObject    handle to runInterp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MTEX = getappdata(handles.output,'MTEX');
filterVal = str2double(get(handles.filterVal,'string'));
[ebsd,grains,DIC] = straininterp7_7GUI(MTEX.DIC.ProjFileName,MTEX.DIC.ProjFilePathName,MTEX.EBSD.ProjFileName,MTEX.EBSD.ProjFilePathName,MTEX.Cal.ProjFileName,MTEX.Cal.ProjFilePathName,filterVal,handles);
% setappdata(handles.axes1,'ebsd',ebsd)
% setappdata(handles.axes1,'grains',grains)
% setappdata(handles.axes1,'DIC',DIC)

% --- Executes on button press in browseTXT.
function browseTXT_Callback(hObject, eventdata, handles)
% hObject    handle to browseTXT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[DIC.ProjFileName,DIC.ProjFilePathName] = uigetfile('*.txt','read text file','MultiSelect', 'on');
setappdata(handles.output,'DIC',DIC)
%set(handles.DICroot,'String',DIC.ProjFileName(1))
set(handles.DICroot,'String',DIC.ProjFilePathName)


function DICroot_Callback(hObject, eventdata, handles)
% hObject    handle to DICroot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DICroot as text
%        str2double(get(hObject,'String')) returns contents of DICroot as a double


% --- Executes during object creation, after setting all properties.
function DICroot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DICroot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in runStrain.
function runStrain_Callback(hObject, eventdata, handles)
% hObject    handle to runStrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DIC = getappdata(handles.output,'DIC');

DIC.window = round(str2double(get(handles.winSize,'string')));
DIC.step = round(str2double(get(handles.stepSize,'string')));


DIC.SW(1) = round(str2double(get(handles.SWx,'string')));
DIC.SW(2) = round(str2double(get(handles.SWy,'string')));
DIC.SWpix = DIC.window + DIC.SW*DIC.step;%

DIC = straincalc6_3GUI(DIC,handles);
%setappdata(handles.output,'DIC',DIC)


function winSize_Callback(hObject, eventdata, handles)
% hObject    handle to winSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of winSize as text
%        str2double(get(hObject,'String')) returns contents of winSize as a double


% --- Executes during object creation, after setting all properties.
function winSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to winSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stepSize_Callback(hObject, eventdata, handles)
% hObject    handle to stepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stepSize as text
%        str2double(get(hObject,'String')) returns contents of stepSize as a double


% --- Executes during object creation, after setting all properties.
function stepSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SWx_Callback(hObject, eventdata, handles)
% hObject    handle to SWx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SWx as text
%        str2double(get(hObject,'String')) returns contents of SWx as a double
if get(handles.linkxy,'Value') ==1
    SW = str2double(get(handles.SWx,'string'));
    set(handles.SWy,'string',SW)
end

% --- Executes during object creation, after setting all properties.
function SWx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SWx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SWy_Callback(hObject, eventdata, handles)
% hObject    handle to SWy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SWy as text
%        str2double(get(hObject,'String')) returns contents of SWy as a double
if get(handles.linkxy,'Value') ==1
    SW = str2double(get(handles.SWy,'string'));
    set(handles.SWx,'string',SW)
end

% --- Executes during object creation, after setting all properties.
function SWy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SWy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in overlapSW.
function overlapSW_Callback(hObject, eventdata, handles)
% hObject    handle to overlapSW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of overlapSW
overlapVal = get(handles.overlapSW,'Value');
if overlapVal == 1
    set(handles.overlapSW,'Value',0)
else
    set(handles.overlapSW,'Value',1)
end


% --- Executes during object creation, after setting all properties.
function overlapSW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to overlapSW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
overlapVal = get(handles.overlapSW,'Value');
if overlapVal == 1
    set(handles.overlapSW,'Value',0)
else
    set(handles.overlapSW,'Value',1)
end


% --- Executes on button press in linkxy.
function linkxy_Callback(hObject, eventdata, handles)
% hObject    handle to linkxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of linkxy
if get(handles.linkxy,'Value')==1
    SW = str2double(get(handles.SWx,'string'));
    set(handles.SWy,'string',SW)
end
