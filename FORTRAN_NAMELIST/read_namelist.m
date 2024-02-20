function data=read_namelist(filename,namelist)
%data=READ_NAMELIST(filename,namelist) Returns values from FORTRAN namelist
%   This function reads a FORTRAN input namelist and returns the values as
%   the fields of a structure.  Multidimensional arrays have their indicies
%   scaled and shifted to fit the matlab numbering scheme.  Alternatively
%   the user may specify an open file id.
%
%   Note:  At this time the function does not support the collon modifier
%   and variables such as TEST(1:,5) = 1 2 3 4 will be overlooked.
%
%   Example:
%       data=read_namelist('input.test','INDATA');
%       fid = fopen('input.test','r');
%       data=read_namelist(fid,'INDATA');
%       fclose(fid);
%
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.5
%   Date:           9/13/12

% Open a text file
if ischar(filename)
    fid=fopen(filename,'r');
else
    fid = filename;
end
% Find the namelist section
line=fgetl(fid);
while ~feof(fid) && isempty(strfind(line,namelist))
    line=fgetl(fid);
end
if feof(fid)
    disp(['ERROR: Namelist: ' namelist ' not found!']);
    data=-1;
    return
end
data.temp=-1;
% Now read the namelist
line=fgetl(fid);
total_line=[];
while (~strncmpi(strtrim(line),'/',1) && (~strncmpi(strtrim(upper(line)),'&END',4) && ~strncmpi(strtrim(upper(line)),'$END',4)))
    % Get rid of commentary lines
    line = strtrim(line);
    dex = strfind(line,'!');
    if ~isempty(dex)
        if (dex(1) == 1)
            line = '';
        else
            line = line(1:dex(1)-1);
        end
    end
    if length(line) > 1
        total_line=[total_line '  ' line];
    end
    % Get next nonblank line.
    line=fgetl(fid);
    while isempty(line)
        line=fgetl(fid);
        strtrim(upper(line));
    end
end
% Close the text file
if ischar(filename)
    fclose(fid);
end
% Now we extract every declaration
expr=['(([a-zA-Z0-9_]+)|([a-zA-Z0-9_]+\([\+\-0-9 ,:]+\)))\s+='...
    '(((\s+(([0-9]+[*])?(\+|\-)?[0-9.]+E(\+|\-)?[0-9]+)+)+)|'...
    '(\s+\''[a-zA-Z0-9_/.]+\'')|(\s+\''.+\'')'...
    '(\s+[T|F|\.true\.|\.false\.])|'...
    '((\s+([0-9]+[*])?(\+|\-)?[0-9.]+)+))'];  %<- Needs to be modified for comma delimiters
total_line=regexprep(total_line,'([\n|\r])','','ignorecase'); %Get Rid of CR
total_line=regexprep(total_line,' T ',' 1 ','ignorecase'); %Get Rid of T
total_line=regexprep(total_line,' F ',' 0 ','ignorecase'); %Get Rid of F
total_line=regexprep(total_line,' \.true\. ',' 1 ','ignorecase'); %Get Rid of .true.
total_line=regexprep(total_line,' \.false\. ',' 0 ','ignorecase'); %Get Rid of .false.
total_line=regexprep(total_line,',',' ','ignorecase'); % Get Rid of commas (B. Bowler)
% Fix spacing issue
dex = find(total_line == '=');
for i=1:length(dex)
    tstr='=';
    if (total_line(dex(i)-1) ~= ' ')
        tstr=' =';
    end
    if (total_line(dex(i)+1) ~= ' ')
        tstr=[tstr ' '];
    end
    total_line = [total_line(1:dex(i)-1) tstr total_line(dex(i)+1:length(total_line))];
    dex = find(total_line == '=');
end
lines=regexpi(total_line,expr,'match');
% Now we parse each declaration
nlines=numel(lines);
for i=1:nlines
    line=handlestars(lines{i});
    eqdex=strfind(line,'=');
    par1dex=strfind(line,'(');
    par2dex=strfind(line,')');
    comadex=strfind(line,',');
    if isempty(comadex)
        if isempty(par1dex)
            name=lower(strtrim(sscanf(line(1:eqdex-1),'%s')));
            vals=sscanf(line(eqdex+1:numel(line)),'%g')';
            if isempty(vals)
                vals=strtrim(sscanf(line(eqdex+1:numel(line)),'%s'));
                vals=vals(2:numel(vals)-1); %Need to do this to handle ' in string
            end
            data.(name)=vals;
        else
            name=lower(strtrim(sscanf(line(1:par1dex-1),'%s')));
            index=sscanf(line(par1dex+1:par2dex-1),'%g');
            vals=sscanf(line(eqdex+1:numel(line)),'%g')';
            if ~isfield(data,name)
                data.(name)=vals;
                if index >= 1
                    data.([name '_nmlindex'])=index;
                end
            else
                data.(name)=[data.(name) vals];
                data.([name '_nmlindex'])=[data.([name '_nmlindex']) index];
            end
        end
    else
        index_string='%g';
        for j=1:numel(comadex)
            index_string=[index_string ',%g'];
        end
        name=lower(strtrim(sscanf(line(1:par1dex-1),'%s')));
        index=sscanf(line(par1dex+1:par2dex-1),index_string);
        vals=sscanf(line(eqdex+1:numel(line)),'%g')';
        if ~isfield(data,name)
            data.(name)=vals;
            for j=1:numel(comadex)+1
                data.([name '_nmlindex' num2str(j)])=index(j);
                data.([name '_maxorder'])=numel(comadex)+1;
            end
        else
            data.(name)=[data.(name) vals];
            for j=1:numel(comadex)+1
                data.([name '_nmlindex' num2str(j)])=[data.([name '_nmlindex' num2str(j)]) index(j)];
            end
        end
    end
end
% Now we reformulate the arrays
data=rmfield(data,'temp');
names=fieldnames(data);
for i=1:numel(names)
    if isfield(data,[names{i} '_nmlindex'])
        for j=1:numel(data.(names{i}))
            temp(data.([names{i} '_nmlindex'])(j))=data.(names{i})(j);
        end
        data.(names{i})=temp;
        data=rmfield(data,[names{i} '_nmlindex']);
    elseif isfield(data,[names{i} '_nmlindex1'])
        test=0;
        for j=1:data.([names{i} '_maxorder']);
            if isfield(data,[names{i} '_nmlindex' num2str(j)])
                test=1;
            else
                test=0;
            end
        end
        if test
            minmax=zeros(2,data.([names{i} '_maxorder']));
            for j=1:data.([names{i} '_maxorder'])
                minmax(1,j)=min(data.([names{i} '_nmlindex' num2str(j)]));
                minmax(2,j)=max(data.([names{i} '_nmlindex' num2str(j)]));
            end
            arraysize=minmax(2,1)-minmax(1,1)+1;
            for j=2:data.([names{i} '_maxorder'])
                arraysize=[arraysize minmax(2,j)-minmax(1,j)+1];
            end
            temp=zeros(arraysize);
            % Now we redo the _nmlindex arrays to get the proper ordering
            for j=1:data.([names{i} '_maxorder'])
                if minmax(1,j)< 1
                    data.([names{i} '_nmlindex' num2str(j)])=...
                        data.([names{i} '_nmlindex' num2str(j)])-minmax(1,j)+1;
                end
            end
            % Create the new array by creating a string to index multiple
            % indexes
            for j=1:numel(data.(names{i}))
                % Create a executable string
                exestring=['temp('];
                for k=1:data.([names{i} '_maxorder'])
                    exestring=[exestring 'data.' names{i} '_nmlindex' num2str(k)...
                        '(' num2str(j) '),'];
                end
                exestring=[exestring(1:numel(exestring)-1) ')=data.' names{i}...
                    '(' num2str(j) ');'];
                eval(exestring);
            end
            data.(names{i})=temp;
            % Now cleanup fields
            for k=1:data.([names{i} '_maxorder'])
                data=rmfield(data,[names{i} '_nmlindex' num2str(k)]);
            end
            data=rmfield(data,[names{i} '_maxorder']);
        end
    end
end
return
end

function output=handlestars(input)
% This function handles the multiple field references made in an input
% namelist such as 128*5.5.
output=input;
stardex=strfind(output,'*');
while ~isempty(stardex)
        spdex1=strfind(output(1:stardex(1)),' ');
        spdex1=spdex1(numel(spdex1));
        numvals=sscanf(output(spdex1:stardex(1)-1),'%g');
        spdex2=strfind(output(stardex(1):numel(output)),' ');
        front=output(1:spdex1);
        if isempty(spdex2)
            spdex2=numel(output);
            val=sscanf(output(stardex(1)+1:numel(output)),'%g');
            back='';
        else
            spdex2=spdex2(1)+stardex(1);
            val=sscanf(output(stardex(1)+1:spdex2-1),'%g');
            back=output(spdex2:numel(output));
        end
        output=front;
        for i=1:numvals
            output=[output ' ' num2str(val)];
        end
        output=[output ' ' back];
        stardex=strfind(output,'*');
end
return
end