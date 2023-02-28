function makeStimMovie(stimMat, fileName)

maxVal = 15;
close all

nFrames = size(stimMat, 3);

% Preallocate movie structure.
writerObj = VideoWriter(fileName);

open(writerObj);

set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

% Create movie.

for ii = 1:nFrames 
   plotMidFrameG4(fliplr(stimMat(:,:,ii)), maxVal)
   frame = getframe;
   writeVideo(writerObj,frame);
end

% Create AVI file.
close(writerObj);