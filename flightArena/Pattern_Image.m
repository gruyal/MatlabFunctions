function Pattern_Image(Pats,gs_val,filename,format)

%parameters
savedir = 'C:\matlabroot\Pattern Images\';
color = [0 1 0]; %green
magn = 2; %magnification
default_format = 'png'; %png, gif, or avi
skip_first_frame = 1; %1 to skip the first frame of the pattern
rot180 = 1; %1 if arena is mounted upside-down

if ~exist(savedir,'dir')
    mkdir(savedir);
end

if nargin<4
    [filename, filedir] = uigetfile('Pattern*.mat','Select pattern file');
    load([filedir filename]);
    filename = filename(1:end-4); %remove format from string
    Pats = pattern.Pats(:,:,1+skip_first_frame:end);
    if rot180==1
        Pats = rot90(Pats,2);
    end
    gs_val = pattern.gs_val;
    format = default_format;
end
file = [savedir filename '.' format];

%increase size of gif and convert to 0-1 scale
S = zeros(size(Pats).*[magn magn 1]);
for f = 1:size(Pats,3)
    S(:,:,f) = kron(Pats(:,:,f)/(2^gs_val-1),ones(magn)); 
end

%convert gif to desired color
S = permute(S,[1 2 4 3]);
color = permute(color,[1 3 2]);
color = repmat(color,size(S));
S = repmat(S,[1 1 3 1]);
S = S.*color; 

switch format
    case 'gif'
        for f = 1:size(S,4)
          [frame,map] = rgb2ind(S(:,:,:,f),2^gs_val,'nodither');
          if f == 1
            imwrite(frame,map,file,'Loop',Inf,'Delay',0);
          else
            imwrite(frame,map,file,'WriteMode','append','Delay',0);
          end
        end

    case 'avi'
        num_loops = 3;
        outputVideo = VideoWriter(file);
        outputVideo.FrameRate = 10;
        open(outputVideo)
        for f = repmat(1:size(S,4),1,num_loops)
           writeVideo(outputVideo,S(:,:,:,f))
        end
        close(outputVideo)
        
    case 'png'
        [frame,map] = rgb2ind(S(:,:,:,1),2^gs_val,'nodither');
        imwrite(frame,map,file);
end

disp(['pattern saved as "' file '"'])

end