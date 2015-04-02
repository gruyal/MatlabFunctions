function maskIm = generateBaseMask(maskType, maskSize)


% maskIm = generateBaseMask(maskType, maskSize)
%
% This function generates a basic mask to be used with the grating patterns
% that are generated by generateGratingFrame. The mask can be either a
% circle or a square (maskType, given as a string), and will be set in the
% middle of the 99X99 matrix with a diameter/vertex of 2XmaskSize+1 pixels
%
% Note! due to rotation symmetry circles are not available in all sizes
%
% OUTPUT
%
% maskIm -      225X225 matrix with the mask centered <(96+16)X2 + 1 full
% arena + biggest maskSize>


% these reduce or eliminate asymmetries in 45 degree rotations (radius,
% linspacesteps)
circleComb = [2, 11; 3, 15; 4, 41; 5, 29; 7, 45; 9, 77; 10, 131; 12, 101; 15, 133; 17, 141]; 

matSiz = 225;
cen=ceil(matSiz/2); 
maskIm = zeros(cen*2-1);

assert(maskSize > 0, 'maskSize must be positive')
assert(maskSize < 18,'MaskSize cannot exceed 17')


if strcmpi(maskType, 'square')
    
    maskIm(cen-maskSize:cen+maskSize, cen-maskSize:cen+maskSize)=1;
    
elseif strcmpi(maskType, 'circle')
    
    combInd = find(circleComb(:, 1)== maskSize, maskSize);
    
    if isempty(combInd)
        fprintf('relevant sizes are: %d, %d, %d, %d, %d, %d, %d, %d, %d\n', circleComb(:,1)) 
        error('Cannot generate a rotation symmetric circle at this size')
        
    end
    
    base = linspace(0, 2*pi, circleComb(combInd, 2));
    radius = maskSize;

    rx = round(sin(base)*radius + cen);
    ry = round(cos(base)*radius + cen);
    
    rx = [rx'; rx(1)];
    ry = [ry'; ry(1)];
    
    
%    maskIm = poly2mask(rx,ry, matSiz, matSiz);
    
    

    for ii=1:length(rx)
        maskIm(rx(ii),ry(ii)) =1;
    end

    [rr, cc] = find(maskIm);
    urr = unique(rr);

    for ii=1:length(urr)
        tempind = find(rr == urr(ii));
        maskIm(urr(ii), min(cc(tempind)):max(cc(tempind))) = 1;    
    end

else
    error('mask type must be either square or circle')

end





end
