
function [ob index]=convertTDMS(filename)

%convertTDMS - function to convert Labview TDMS data files into .mat files.
%   If called with no input, user selects from file open dialog box.  A
%   .mat file with the same filename as the TDMS file is automatically
%   written to the same directory (warning - will over write .mat file of 
%   the same name. 
%
%   TDMS format is based on information provided by National Instruments
%   at:    http://zone.ni.com/devzone/cda/tut/p/id/5696
%
% convertTDMS(filename)
%
%       Inputs:
%               filename - filename to be converted.  If not supplied, user
%                 is provided dialog box to open file.  Can be a cell array
%                 of files for bulk conversion.
%
%       Outputs:
%               ob - structure with all of the data objects
%               index - index of the information in ob
%

%---------------------------------------------
%Brad Humphreys - v1.0 2008-04-23
%ZIN Technologies
%---------------------------------------------

conver='1.0';    %conversion version number

startingdir=cd; %Get the starting directory

if nargin==0                %If file is not provided,prompt user
    [filename,pathname,filterindex] = uigetfile({'*.*',  'All Files (*.*)'},'Choose File');
    filename=fullfile(pathname,filename);
end

%Create filelist
if iscell(filename)         %If it is a group of files
    numfiles=max(size(filename));
    infilename=filename;
else
    numfiles=1;             %If only one file name is provided
    infilename{1}=filename;
end



for fnum=1:numfiles
    disp([datestr(now,13) ' Begining conversion of:  '  infilename{fnum}])

    bytesInfile = dir([infilename{fnum}, '.tdms']);    %Get size of file.  Needed for later estimation of variable size.
    filesize = bytesInfile.bytes;

    %Initialize variables for each file conversion
    index.names=[];
    ob=[];
    firstRawData=1;
    segCnt=0;

    fid=fopen([infilename{fnum}, '.tdms']);
    
    while ~feof(fid)
        
        %Look for TDM tag
        Ttag=fread(fid,1,'uint8');
        Dtag=fread(fid,1,'uint8');
        Stag=fread(fid,1,'uint8');
        mtag=fread(fid,1,'uint8');
        if Ttag==84 && Dtag==68 && Stag==83 && mtag==109
            segCnt=segCnt+1;     
        else
            if ~(isempty(Ttag) && isempty(Dtag) && isempty(Stag) && isempty(mtag))  %On last read, all will be empty
                error('Unexpected bit stream in file: %s  at segment %d (dec: %d %d %d %d)',infilename{fnum}, segCnt+1, Ttag, Dtag, Stag, mtag);
            end
        end

        %ToC Field
        Toc=fread(fid,1,'uint32');
        kTocMetaData=bitget(Toc,2);
        kTocNewObject=bitget(Toc,3);
        kTocRawData=bitget(Toc,4);

        %Segment
        vernum=fread(fid,1,'uint32');
        segLength=fread(fid,1,'uint64');
        metaRawLength=fread(fid,1,'uint64');

 %% Process Meta Data
        if kTocMetaData                                     %If there is meat data in this segment
            numNewObjInSeg=fread(fid,1,'uint32');
            for q=1:numNewObjInSeg

                obLength=fread(fid,1,'uint32');             %Get the length of the objects name
                obname=[fread(fid,obLength,'uint8=>char')]';%Get the objects name

                %Fix Object Name
                if strcmp(obname,'/')
                    obname='Root';
                else
                    obname=fixcharformatlab(obname);
                end

                if ~isfield(ob,obname)                         %If the object does not already exist, create it
                    index.names{end+1}=obname;
                    ob.(obname)=[];
                    obnum=max(size(index.names));               %Get the object number
                else                                            %if it does exist, get it's index and object number
                    obnum=find(strcmp(index.names,obname)==1,1,'last');
                end

                rawdataindex=fread(fid,1,'uint32');             %Get the raw data Index
                if rawdataindex==0                              %No raw data assigned to this object in this segment
                    index.entry(obnum)=0;
                    index.dataType(obnum)=1;
                    index.arrayDim(obnum)=0;
                    index.nValues(obnum)=0;
                    index.byteSize(obnum)=0;
                elseif rawdataindex+1==2^32                     %Objects raw data index matches previous index - no changes
                else                                            %Get new object information
                    index.entry(obnum)=rawdataindex;
                    index.dataType(obnum)=fread(fid,1,'uint32');
                    index.arrayDim(obnum)=fread(fid,1,'uint32');
                    index.nValues(obnum)=fread(fid,1,'uint64');
                    if index.dataType(obnum)==32                %If the datatype is a string
                        index.byteSize(obnum)=fread(fid,1,'uint64');
                    else
                        index.byteSize(obnum)=0;
                    end
                end


                numProps=fread(fid,1,'uint32');
                for p=1:numProps
                    propNameLength=fread(fid,1,'uint32');
                    propsName=[fread(fid,propNameLength,'uint8=>char')]';
                    propsDataType=fread(fid,1,'uint32');
                    propExists=isfield(ob.(obname),propsName);
                    dataExists=isfield(ob.(obname),'data');

                    if dataExists                                               %Get number of data samples for the object in this segment
                        %nsamps=max(size(ob.(obname).data));
                        
                        
                        nsamps=ob.(cname).nsamples+1;
                        
                    else
                        nsamps=0;
                        estNumSeg=floor(filesize/(20+segLength)*1.2);            %Estimate # of Segements.  20 is the number of bytes prior to the segLength Read
                    end

                    if propsDataType==32                                         %If the data type is a string
                        propsValueLength=fread(fid,1,'uint32');
                        propsValue=fread(fid,propsValueLength,'uint8=>char');
                        if propExists
                            cnt=ob.(obname).(propsName).cnt+1;
                            ob.(obname).(propsName).cnt=cnt;
                            ob.(obname).(propsName).value{cnt}=propsValue;
                            ob.(obname).(propsName).samples(cnt)=nsamps;
                        else
                            ob.(obname).(propsName).cnt=1;
                            ob.(obname).(propsName).value{estNumSeg}=NaN;
                            ob.(obname).(propsName).samples(estNumSeg)=0;
                            ob.(obname).(propsName).value{1}=propsValue;
                            ob.(obname).(propsName).samples(1)=nsamps;
                        end
                    else                                                        %Numeric Data type
                        if propsDataType==68                                     %If the data type is a timestamp
                            tsec=fread(fid,1,'uint64')/2^64+fread(fid,1,'uint64');   %time since Jan-1-1904 in seconds
                            propsValue=tsec/86400+695422-5/24;                   %/864000 convert to days; +695422 days from Jan-0-0000 to Jan-1-1904
                        else
                            matType=LV2MatlabDataType(propsDataType);
                            propsValue=fread(fid,1,matType);
                        end
                        if propExists
                            cnt=ob.(obname).(propsName).cnt+1;
                            ob.(obname).(propsName).cnt=cnt;
                            ob.(obname).(propsName).value(cnt)=propsValue;
                            ob.(obname).(propsName).samples(cnt)=nsamps;
                        else
                            ob.(obname).(propsName).cnt=1;
                            ob.(obname).(propsName).value(estNumSeg)=NaN;
                            ob.(obname).(propsName).samples(estNumSeg)=0;
                            ob.(obname).(propsName).value(1)=propsValue;
                            ob.(obname).(propsName).samples(1)=nsamps;
                        end
                    end

                end

            end
        end

        %% Process Raw Data
        if kTocRawData                                                          %Process Raw Data
            if firstRawData                                                     %If first raw data, preallocate data arrays
                firstRawData=0;
                estNumSeg=filesize/(20+segLength);                              %Estimate # of Segements.  20 is the number of bytes prior to the segLength Read
                                                                               
                for b=1:max(size(index.names))
                    cname=cell2mat(index.names(b));
                    nEstSamples=floor(1.2*estNumSeg*index.nValues(b));
                    if index.dataType(b)~=32                                    %If the data is numeric type
                        ob.(cname).data(1:(nEstSamples),1)=NaN;
                    else                                                        %If the data is string type
                        ob.(cname).data{nEstSamples}=NaN;
                    end
                    ob.(cname).nsamples=0;
                end
            end

            for r=1:max(size(index.names))                                      %Loop through the index
                nvals=index.nValues(r);
                if nvals>0
                    [data cnt]=fread(fid,nvals,LV2MatlabDataType(index.dataType(r)));
                    cname=cell2mat(index.names(r));
                    ssamples=ob.(cname).nsamples;
                    if index.dataType(r)~=32                                    %If the data is numeric type
                        ob.(cname).data(ssamples+1:ssamples+cnt,1)=data;
                    else                                                        %If the data is string type
                        ob.(cname).data{ssamples+1:ssamples+cnt,1}=data;
                    end
                    ob.(cname).nsamples=ssamples+cnt;
                end
            end
        end
    end


    %% Clean up preallocated arrays   (preallocation required for speed)
    for y=1:max(size(index.names))
        cname=cell2mat(index.names(y));
        nsamples=ob.(cname).nsamples;
        if nsamples>0       %Remove any excess from preallocation of data
            if index.dataType(y)~=32
                ob.(cname).data=ob.(cname).data(1:nsamples,1);                  %If the data is numeric type
            else
                ob.(cname).data=ob.(cname).data{1:nsamples,1};                  %If the data is string type
            end
            
            proplist=fieldnames(ob.(cname));    %Remove any excess from preallocation of properties
            for isaac=1:size(proplist,1)
                if isfield(ob.(cname).(proplist{isaac}),'cnt')
                    cnt=ob.(cname).(proplist{isaac}).cnt;
                    ob.(cname).(proplist{isaac}).value=ob.(cname).(proplist{isaac}).value(1:cnt);
                    ob.(cname).(proplist{isaac}).samples=ob.(cname).(proplist{isaac}).samples(1:cnt);
                    ob.(cname).(proplist{isaac})=rmfield(ob.(cname).(proplist{isaac}),'cnt');
                end
            end
            
        end
    end

    

        %% Create Output filenames
        [pathstr,namestr]=fileparts(infilename{fnum});
        namestr=fixcharformatlab(namestr);
        %savefilename=fullfile(pathstr,namestr);

        if ~isempty(pathstr)
            cd(pathstr)
        end
        fclose(fid);
        
        %Save Output Data
        ob.index=index;
        ob.conver=conver;
        save(namestr,'-struct','ob')
        disp([datestr(now,13) ' Completed conversion of:  '  infilename{fnum}])


    end

    cd(startingdir);

end

    function  fixedtext=fixcharformatlab(textin)
        %Private Function to remove all text that is not MATLAB variable name
        %compatible
        textin=strrep(textin,'''','');
        textin=strrep(textin,'/Untitled/','');
        textin=strrep(textin,'/','.');
        textin=strrep(textin,'-','');
        textin=strrep(textin,'?','');
        textin=strrep(textin,' ','');
        textin=strrep(textin,'.','');
        textin=strrep(textin,'[','_');
        textin=strrep(textin,']','');
        textin=strrep(textin,'%','');
        textin=strrep(textin,'#','');
        textin=strrep(textin,'(','');
        textin=strrep(textin,')','');
        fixedtext=textin;

    end

    function matType=LV2MatlabDataType(LVType)
        %Cross Refernce Labview TDMS Data type to MATLAB

        switch LVType
            case 1   %tdsTypeVoid
                matType='';
            case 2   %tdsTypeI8
                matType='int8';
            case 3   %tdsTypeI16
                matType='int32';
            case 4   %tdsTypeI32
                matType='int32';
            case 5   %tdsTypeI64
                matType='int64';
            case 6   %tdsTypeU8
                matType='uint8';
            case 7   %tdsTypeU16
                matType='uint16';
            case 8   %tdsTypeU32
                matType='uint32';
            case 9   %tdsTypeU64
                matType='uint64';
            case 10  %tdsTypeSingleFloat
                matType='float64';
            case 11  %tdsTypeDoubleFloat
                matType='float64';
            case 12  %tdsTypeExtendedFloat
                matType='';
            case 32  %tdsTypeString
                matType='char';
            case 33  %tdsTypeBoolean
                matType='bit1';
            case 68  %tdsTypeTimeStamp
                matType='bit224';
            otherwise
                matType='';
        end

    end

