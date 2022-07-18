function DICEebsd
%% Define input file

%get file name and path
[ctffile,path] = uigetfile('*.crc');
% [ctffile,path] = uigetfile({'*.cpr;*.crc'});

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [3.6 3.6 3.6], 'mineral', 'Ni-superalloy', 'color', [0.53 0.81 0.98])};

% plotting convention
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','intoPlane');

%% Import the Data

% create an EBSD variable containing the data
ebsd = EBSD.load([path filesep ctffile],CS,'interface','crc',...
  'convertEuler2SpatialReferenceFrame');

%% Sort out input for DICE

%remove non-indexed points
ebsd = ebsd('indexed');

%reconstruct grains
[grains,ebsd.grainId] = calcGrains(ebsd,'angle',15*degree);

%smooth grains
grains = smooth(grains,5);
% gB = grains.boundary

save([path filesep ctffile(1:end-4)],'ebsd','grains')
