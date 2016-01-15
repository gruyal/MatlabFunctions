function maskIm = makeRectangleMask(maskhW, maskhH, matSize)

% function maskIm = makeSquareMask(maskhE, matSize)
% this function generates a square mask to be used by generateBaseMask
%
% INPUTS
% maskhW/H - half width/height for specific mask (rectangle will have an width/height of
% 2*maskhW/H+1. Width is along the arenas X dimension and height is arena Y
% matSize - Size of the matrix that the mask will be placed in its center
%
% NOTE! matSize should be odd so that rotations of the image will do
% minimal distortion. 
%
% OUTPUT
% maskIm - matSize X matSize matrix with a maskE edge of 1's in the middle


assert(maskhW > 0, 'mask width must be a positive number')
assert(maskhH > 0, 'mask height be a positive number')
assert(matSize > 0,'matSize should be a positive number')
assert(matSize/2 ~= floor(matSize/2), 'matSize should be odd to avoid rotation distortions')

cen=ceil(matSize/2); 
maskIm = zeros(matSize);


maskIm(cen-maskhH:cen+maskhH, cen-maskhW:cen+maskhW) = 1;



end