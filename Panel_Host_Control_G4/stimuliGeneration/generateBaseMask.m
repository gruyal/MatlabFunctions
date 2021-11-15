function maskIm = generateBaseMask(maskStruct)


% maskIm = generateBaseMask(maskStruct)
%
% This function generates a basic mask to be used with the grating patterns
% that are generated by generateGratingFrame. The mask can be either a
% circle, annulus, square, or rectangle (maskType, given as a string), and will be set in the
% middle of the matSizXmatSiz matrix with a diameter/Edge of 2XmaskSize+1 pixels
%
% Note! for annulus maskSize should be a 2 element vector. Bigger number
% would be used for outer and smaller for inner
% Note! due to rotation symmetry circles are not available in all sizes
%
% INPUT
% maskStruct -      strcture that contains the following fileds:
% .type             mask type (string). options are 'circle', 'square',
%                   'annulus', and 'rectangle'
% .radius -         radius. depending on type this is either a single
%                   number or a 1X2 vector
% .ori -            orientation. Only applied to rectangle. allows to use
%                   generateBarFrameByInds by circomventing the imrotate use within
%                   createProtocol. 0-7 multiples of 45 degrees
%
% OUTPUT
%
% maskIm -  445X445 matrix with the mask centered <(96+16)X2 + 1 full
% arena + biggest maskSize>


% these reduce or eliminate asymmetries in 45 degree rotations (radius,
% linspacesteps)

matSiz = 445;

assert(isfield(maskStruct, 'type'), 'maskStruct is missing type')
assert(isfield(maskStruct, 'radius'), 'maskStruct is missing radius')
maskSize = maskStruct.radius;
maskType = maskStruct.type;
maskType = lower(maskType);

switch maskType

    case 'square'

        maskIm = makeSquareMask(maskSize, matSiz);

    case 'circle'

        maskIm = makeCircleMask(maskSize, matSiz);

    case 'annulus'

        assert(length(maskSize) == 2, 'When type is annulus - maskSize should be a 2 element vector (R outer, R inner)')
        maskIm = makeAnnulusMask(max(maskSize), min(maskSize), matSiz);

    case 'rectangle'
        assert(length(maskSize) == 2, 'When type is rectangle - maskSize should be a 2 element vector (width, height)')
        if isfield(maskStruct, 'ori')
            ori = maskStruct.ori;
            maskIm = makeRectangleMask(maskSize(1), maskSize(2), matSiz, ori);
        else
            maskIm = makeRectangleMask(maskSize(1), maskSize(2), matSiz);
        end

    otherwise

        error('mask type must be either square, circle, annulus, or rectangle')

end



end