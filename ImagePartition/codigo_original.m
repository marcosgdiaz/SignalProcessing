% Demo to use normxcorr2 to find a template (a white onion)
% in a larger image (of a pile of vegetables)
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
clear;  % Erase all existing variables.
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

% Read in a standard MATLAB color demo image.
folder = '/';
baseFileName = 'test01i.jpg';
% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
if ~exist(fullFileName, 'file')
	% Didn't find it there.  Check the search path for it.
	fullFileName = baseFileName; % No path this time.
	if ~exist(fullFileName, 'file')
		% Still didn't find it.  Alert user.
		errorMessage = sprintf('Error: %s\\%s does not exist.', folder, fullFileName);
		uiwait(warndlg(errorMessage));
		return;
	end
end
rgbImage = imread(fullFileName);
% Get the dimensions of the image.  numberOfColorBands should be = 3.
[rows, columns, numberOfColorBands] = size(rgbImage)
% Display the original image.
subplot(2, 3, 1);
imshow(rgbImage, []);
axis on;
caption = sprintf('Original Color Image, %d rows by %d columns.', rows, columns);
title(caption, 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0, 0, 1, 1]);

% Convert to binary
if numberOfColorBands > 1
	binaryImage = rgbImage(:,:,1) > 200;
else
	binaryImage = rgbImage > 200;
end
% Display the original image.
subplot(2, 3, 2);
imshow(binaryImage, []);
axis on;
caption = sprintf('Original Binary Image, %d rows by %d columns.', rows, columns);
title(caption, 'FontSize', fontSize);

% Clean up a bit.
% binaryImage = imclearborder(binaryImage);
binaryImage = imfill(~binaryImage, 'holes');
% Display the image.
subplot(2, 3, 3);
imshow(binaryImage, []);
axis on;
caption = sprintf('Cleaned Binary Image, %d rows by %d columns.', rows, columns);
title(caption, 'FontSize', fontSize);
drawnow; % Force it to refresh.

labeledImage = bwlabel(binaryImage, 8);     % Label each blob so we can make measurements of it
coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels

subplot(2, 3, 4);
imshow(labeledImage, []);
title('Labeled Image, from bwlabel()', 'FontSize', fontSize);
subplot(2, 3, 5);
imshow(coloredLabels);
caption = sprintf('Pseudo colored labels, from label2rgb().\nBlobs are numbered from top to bottom, then from left to right.');
title(caption, 'FontSize', fontSize);
drawnow; % Force it to refresh.

% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
blobMeasurements = regionprops(labeledImage, 'all');
numberOfBlobs = size(blobMeasurements, 1);

% Now I'll demonstrate how to select certain blobs based using the ismember function.
% Let's say that we wanted to find only those blobs
% with an intensity between 150 and 220 and an area less than 2000 pixels.
% This would give us the three brightest dimes (the smaller coin type).
%allPerimeters = [blobMeasurements.Perimeter]
allAreas = [blobMeasurements.Area];
%allCircularities = allPerimeters .^ 2 ./ (4 * pi * allAreas);
%sortedCirc = sort(allCircularities, 'Ascend')
	
% Get a list of the blobs that meet our criteria and we need to keep.
%allowablecircIndexes = (allCircularities < 6);
allowableAreaIndexes = allAreas > 100; % Take the big objects.
keeperIndexes = find(allowableAreaIndexes);
% Extract only those blobs that meet our criteria, and
% eliminate those blobs that don't meet our criteria.
% Note how we use ismember() to do this.
keeperBlobsImage = ismember(labeledImage, keeperIndexes);
% Re-label with only the keeper blobs kept.
labeledImage = bwlabel(keeperBlobsImage, 8);     % Label each blob so we can make measurements of it
% Now we're done.  We have a labeled image of blobs that meet our specified criteria.
subplot(2, 3, 6);
imshow(labeledImage, []);
title('"Keeper" blobs', 'FontSize', fontSize);



