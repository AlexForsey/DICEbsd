% This script reads in an image from a Davis text file and loads it into
% the DIC structure in DICE

%load in image .txt file from Davis. Need to specify filename and location
data = importdata('/Users/Alexander/Downloads/B00001.txt');

% load in DICE project file. Need to specify filename and location
load('/Users/Alexander/Downloads/Strain_SW_1_1_B00101.mat')

clear xStep yStep x0 y0 x y spaces

%find spaces in header
spaces = find(isspace(data.textdata{:}));

%data for x coords
xStep = str2double(data.textdata{:}(spaces(6)+1:spaces(7)-1));
x0 = str2double(data.textdata{:}(spaces(7)+1:spaces(8)-1));

%data for y coords
yStep = str2double(data.textdata{:}(spaces(10)+1:spaces(11)-1));
y0 = str2double(data.textdata{:}(spaces(11)+1:spaces(12)-1));

%build coordinates for image
imSize = size(data.data);
x = x0:xStep:(x0+(xStep*imSize(2)));
y = y0:yStep:(y0+(yStep*imSize(1)));

%plot scaled image
figure(5)
imagesc(x,y,data.data)
axis image

DIC.image.I = data.data;
DIC.image.x = x;
DIC.image.y = y;

%save DICE project file
save('Strain_SW_1_1_B00101_image.mat','DIC')