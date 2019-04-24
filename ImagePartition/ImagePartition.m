function [] = ImagePartition(FullFileName, Umbral, Area, Area2 )
clc;    % Clear the command window.
 % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
 % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 14;

% Check that user has the Image Processing Toolbox installed.
hasIPT = license('test', 'image_toolbox');
if ~hasIPT
	% User does not have the toolbox installed.
	message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
	reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
	if strcmpi(reply, 'No')
		% User said No, so exit.
		return;
	end
end

%Get the full filename, with path prepended.
if ~exist(FullFileName, 'file')
	% Didn't find it there.  Check the search path for it.
	errorMessage = sprintf('Error: %s does not exist.', FullFileName);
	uiwait(warndlg(errorMessage));
	return;
end

rgbImage = imread(FullFileName);
[rows, columns, numberOfColorBands] = size(rgbImage);
% Convert to binary
if numberOfColorBands > 1
    grayImage = rgb2gray(rgbImage);
	binaryImage = grayImage > Umbral;
else
	binaryImage = rgbImage > Umbral;
end

subplot(1,4,1)
imshow(rgbImage,[])
% axis on;
caption=sprintf('Imagen original');
title(caption,'Fontsize',fontSize);
set(gcf,'units','normalized ','OuterPosition',[0 0 1 1])
drawnow;

subplot(1,4,2)
imshow(binaryImage,[])
caption=sprintf('Imagen binaria');
title(caption,'Fontsize',fontSize);
drawnow;

% Clean up a bit.
%binaryImage = imclearborder(binaryImage);
binaryImage = imfill(binaryImage, 'holes');
subplot(1,4,3)
imshow(binaryImage,[])
caption=sprintf('Imagen binaria rellenada');
title(caption,'Fontsize',fontSize);
drawnow;

%Se utilizan funciones de Matlab para quedarnos con los objetas con un area
%mayor que un umbral predefinido
labeledImage = bwlabel(binaryImage, 8);     % Label each blob so we can make measurements of it
coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels
blobMeasurements = regionprops(labeledImage, 'all');
allAreas = [blobMeasurements.Area]
allowableAreaIndexes = allAreas > Area & allAreas < Area2; 
keeperIndexes = find(allowableAreaIndexes);
keeperBlobsImage = ismember(labeledImage, keeperIndexes);
labeledImage = bwlabel(keeperBlobsImage, 8);  
subplot(1, 4, 4);
labeledImage=repmat(uint8(labeledImage),1,1,3).*rgbImage;
labeledImage(:,:,3)=(labeledImage(:,:,3)==0)*140;
imshow(labeledImage, []);
title('Objetos reconocidos', 'FontSize', fontSize);


end

