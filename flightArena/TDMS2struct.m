function data = TDMS2struct(data_folder,data_filename)
% FUNCTION data = TDMS2struct(data_folder, data_filename)

%Specify paths of NI TDM library files
NI_TDM_DLL_Path = 'C:\matlabroot\TDMS_Scripts\MATLAB TDM Example\dev\bin\64-bit\nilibddc.dll';
NI_TDM_H_Path = 'C:\matlabroot\TDMS_Scripts\MATLAB TDM Example\dev\include\64-bit\nilibddc_m.h';
libname = 'nilibddc';

%specify mex compiler
compiler_path = 'C:\TDM-GCC-64';
cc = mex.getCompilerConfigurations;
if exist('cc.Location','var')
    cur_compiler_path = cc.Location;
    if strcmp(cur_compiler_path, compiler_path) == 0
        setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')
        mex -setup
    end
else
    setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')
    mex -setup
end

%Recreate needed property constants defined in nilibddc_m.h
DDC_FILE_NAME					=	'name';
DDC_FILE_DESCRIPTION			=	'description';
DDC_FILE_TITLE					=	'title';
DDC_FILE_AUTHOR					=	'author';
DDC_FILE_DATETIME				=	'datetime';
DDC_CHANNELGROUP_NAME			=	'name';
DDC_CHANNELGROUP_DESCRIPTION	=	'description';
DDC_CHANNEL_NAME				=	'name';

%Check if the paths to 'nilibddc.dll' and 'nilibddc_m.h' have been
%selected. If not, prompt the user to browse to each of the files.
if exist('NI_TDM_DLL_Path','var')==0
    [dllfile,dllfolder]=uigetfile('*dll','Select nilibddc.dll');
    libname=strtok(dllfile,'.');
    NI_TDM_DLL_Path=fullfile(dllfolder,dllfile);
end
if exist('NI_TDM_H_Path','var')==0
    [hfile,hfolder]=uigetfile('*h','Select nilibddc_m.h');
    NI_TDM_H_Path=fullfile(hfolder,hfile);
end

%Set data path. Prompt the user to browse to the path of the TDM or TDMS 
%file to read if not already entered
if nargin==2
    if strcmp(data_filename(end-4:end),'.tdms')==0 || strcmp(data_filename(end-3:end),'.tdm')==0
        data_filename = [data_filename '.tdms'];
    end
    Data_Path = fullfile(data_folder, data_filename);
else
    [filepath,filefolder]=uigetfile({'*.tdm';'*.tdms'},'Select a TDM or TDMS file');
    Data_Path=fullfile(filefolder,filepath);
end

%Load nilibddc.dll (Always call 'unloadlibrary(libname)' after finished using the library)
loadlibrary(NI_TDM_DLL_Path,NI_TDM_H_Path);

%Open the file (Always call 'DDC_CloseFile' when you are finished using a file)
fileIn = 0;
[err,dummyVar,dummyVar,file]=calllib(libname,'DDC_OpenFileEx',Data_Path,'',1,fileIn);

%Read and display file name property
filenamelenIn = 0;
%Get the length of the 'DDC_FILE_NAME' string property
[err,dummyVar,filenamelen]=calllib(libname,'DDC_GetFileStringPropertyLength',file,DDC_FILE_NAME,filenamelenIn);
if err==0 %Only proceed if the property is found
    %Initialize a string to the length of the property value
    pfilename=libpointer('stringPtr',blanks(filenamelen));
    [err,dummyVar,filename]=calllib(libname,'DDC_GetFileProperty',file,DDC_FILE_NAME,pfilename,filenamelen+1);
    setdatatype(filename,'int8Ptr',1,filenamelen);
    data.filename = char(filename.Value);
end

%Read and display file description property
filedesclenIn = 0;
%Get the length of the 'DDC_FILE_DESCRIPTION' string property
[err,dummyVar,filedesclen]=calllib(libname,'DDC_GetFileStringPropertyLength',file,DDC_FILE_DESCRIPTION,filedesclenIn);
if err==0 %Only proceed if the property is found
    %Initialize a string to the length of the property value
    pfiledesc=libpointer('stringPtr',blanks(filedesclen));
    [err,dummyVar,filedesc]=calllib(libname,'DDC_GetFileProperty',file,DDC_FILE_DESCRIPTION,pfiledesc,filedesclen+1);
    setdatatype(filedesc,'int8Ptr',1,filedesclen);
    data.description = char(filedesc.Value);
end

%Read and display file title property
filetitlelenIn = 0;
%Get the length of the 'DDC_FILE_TITLE' string property
[err,dummyVar,filetitlelen]=calllib(libname,'DDC_GetFileStringPropertyLength',file,DDC_FILE_TITLE,filetitlelenIn);
if err==0 %Only proceed if the property is found
    %Initialize a string to the length of the property value
    pfiletitle=libpointer('stringPtr',blanks(filetitlelen));
    [err,dummyVar,filetitle]=calllib(libname,'DDC_GetFileProperty',file,DDC_FILE_TITLE,pfiletitle,filetitlelen+1);
    setdatatype(filetitle,'int8Ptr',1,filetitlelen);
    data.title = char(filetitle.Value);
end

%Read and display file author property
fileauthlenIn = 0;
%Get the length of the 'DDC_FILE_AUTHOR' string property
[err,dummyVar,fileauthlen]=calllib(libname,'DDC_GetFileStringPropertyLength',file,DDC_FILE_AUTHOR,fileauthlenIn);
if err==0 %Only proceed if the property is found
    %Initialize a string to the length of the property value
    pfileauth=libpointer('stringPtr',blanks(fileauthlen));
    [err,dummyVar,fileauth]=calllib(libname,'DDC_GetFileProperty',file,DDC_FILE_AUTHOR,pfileauth,fileauthlen+1);
    setdatatype(fileauth,'int8Ptr',1,fileauthlen);
    data.author = char(fileauth.Value);
end

%Read and display file timestamp property
yearIn = 0;
monthIn = 0;
dayIn = 0;
hourIn = 0;
minuteIn = 0;
secondIn = 0;
msecondIn = 0;
wkdayIn = 0;
[err,dummyVar,year,month,day,hour,minute,second,msecond,wkday]=calllib(libname,'DDC_GetFilePropertyTimestampComponents',file,DDC_FILE_DATETIME,yearIn,monthIn,dayIn,hourIn,minuteIn,secondIn,msecondIn,wkdayIn);
if err==0 %Only proceed if the property is found
    data.timestamp = [num2str(month) '/' num2str(day) '/' num2str(year) ', ' num2str(hour) ':' num2str(minute) ':' num2str(second) ':' num2str(msecond)];
end

%Get channel groups
%Get the number of channel groups
numgrpsIn = 0;
[err,numgrps]=calllib(libname,'DDC_GetNumChannelGroups',file,numgrpsIn);
%Get channel groups only if the number of channel groups is greater than zero
if numgrps>0
	%Initialize an array to hold the desired number of groups
    pgrps=libpointer('int64Ptr',zeros(1,numgrps));
    [err,grps]=calllib(libname,'DDC_GetChannelGroups',file,pgrps,numgrps);
end    

for i=1:numgrps %For each channel group
    %Get channel group name property
    grpnamelenIn = 0;
    [err,dummyVar,grpnamelen]=calllib(libname,'DDC_GetChannelGroupStringPropertyLength',grps(i),DDC_CHANNELGROUP_NAME,grpnamelenIn);
    if err==0 %Only proceed if the property is found
		%Initialize a string to the length of the property value
        pgrpname=libpointer('stringPtr',blanks(grpnamelen));
        [err,dummyVar,grpname]=calllib(libname,'DDC_GetChannelGroupProperty',grps(i),DDC_CHANNELGROUP_NAME,pgrpname,grpnamelen+1);
        setdatatype(grpname,'int8Ptr',1,grpnamelen);
    else
        grpname=libpointer('stringPtr','');
    end
    data.groupname{i} = char(grpname.Value);
        
    %Get channel group description property
    grpdesclenIn = 0;
    [err,dummyVar,grpdesclen]=calllib(libname,'DDC_GetChannelGroupStringPropertyLength',grps(i),DDC_CHANNELGROUP_DESCRIPTION,grpdesclenIn);
    if err==0 %Only proceed if the property is found
		%Initialize a string to the length of the property value
        pgrpdesc=libpointer('stringPtr',blanks(grpdesclen));
        [err,dummyVar,grpdesc]=calllib(libname,'DDC_GetChannelGroupProperty',grps(i),DDC_CHANNELGROUP_DESCRIPTION,pgrpdesc,grpdesclen+1);
        data.groupdesc{i} = char(grpdesc.Value);
    end
        
    
    %Get channels
    numchansIn = 0;
    %Get the number of channels in this channel group
    [err,numchans]=calllib(libname,'DDC_GetNumChannels',grps(i),numchansIn);
    %Get channels only if the number of channels is greater than zero
    if numchans>0
		%Initialize an array to hold the desired number of channels
        pchans=libpointer('int64Ptr',zeros(1,numchans));
        [err,chans]=calllib(libname,'DDC_GetChannels',grps(i),pchans,numchans);
    end
    
    
    for j=1:numchans %For each channel in the channel group
        %Get channel name property
        channamelenIn = 0;
        [err,dummyVar,channamelen]=calllib(libname,'DDC_GetChannelStringPropertyLength',chans(j),DDC_CHANNEL_NAME,channamelenIn);
        if err==0 %Only proceed if the property is found
			%Initialize a string to the length of the property value
            pchanname=libpointer('stringPtr',blanks(channamelen));
            [err,dummyVar,channame]=calllib(libname,'DDC_GetChannelProperty',chans(j),DDC_CHANNEL_NAME,pchanname,channamelen+1);
            setdatatype(channame,'int8Ptr',1,channamelen);
            data.channame{i,j}=char(channame.Value);
        else
            data.channame{i,j}='';
        end
        
        %Get channel data type
        typeIn = 0;
        [err,type]=calllib(libname,'DDC_GetDataType',chans(j),typeIn);
        
        %Get channel values if data type of channel is double (DDC_Double = 10)
        if strcmp(type,'DDC_Double')
            numvalsIn = 0;
            [err,numvals]=calllib(libname,'DDC_GetNumDataValues',chans(j),numvalsIn);
			%Initialize an array to hold the desired number of values
            pvals=libpointer('doublePtr',zeros(1,numvals));
            [err,vals]=calllib(libname,'DDC_GetDataValues',chans(j),0,numvals,pvals);

            setdatatype(vals,'doublePtr',1,numvals);
            
            %Add channel values to a matrix. The comment, #ok<AGROW>, at
            %the end of the line prevents warnings about the matrix needing 
            %to allocate more memory for the added values.
            data.chanvals(i,j,1:length(vals.Value))=(vals.Value); %#ok<AGROW>
        end
    end
end

%Close file
err = calllib(libname,'DDC_CloseFile',file);

%Unload nilibddc.dll
unloadlibrary(libname);

end