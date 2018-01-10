function varargout = VMECedit(varargin)
% VMECEDIT M-file for VMECedit.fig
%      VMECEDIT, is a MATLAB Graphical User Interface (GUI) for creating
%      a VMEC input file (&INPUT namelist).  It allows the user to
%      interactively edit the values in the input file while also seeing
%      visualizations of the various quantities.
%
%      Maintained by: Samuel Lazerson (lazerson@pppl.gov)
%      Version:       0.91

% Edit the above text to modify the response to help VMECedit

% Last Modified by GUIDE v2.5 04-Aug-2010 11:36:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @VMECedit_OpeningFcn, ...
                   'gui_OutputFcn',  @VMECedit_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before VMECedit is made visible.
function VMECedit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VMECedit (see VARARGIN)

% Choose default command line output for VMECedit
handles.output = hObject;
% Create the Data array
handles.data=vmec_namelist_init('indata');
handles.ns_array=[9 19 29 39 49];
handles.ftol_array=[1e-6 1e-8 1e-10 1e-12 1e-15];
handles.am=[0 0 0 0 0 0 0 0 0 0 0];
handles.ai=[0 0 0 0 0 0 0 0 0 0 0];
handles.ac=[0 0 0 0 0 0 0 0 0 0 0];
handles.raxis=[0 0 0 0 0 0 0 0 0 0 0];
handles.zaxis=[0 0 0 0 0 0 0 0 0 0 0];
handles.extcurr=[0 0 0 0];
handles.rbc=zeros(...
    2.*str2double(get(handles.ntor,'String'))+1,...
    str2double(get(handles.mpol,'String'))+1.);
handles.zbs=zeros(...
    2.*str2double(get(handles.ntor,'String'))+1,...
    str2double(get(handles.mpol,'String'))+1.);
handles.nfp_index=5;
handles.ntor_index=6;
handles.mpol_index=9;
% Compose File List
filelist=dir('input.*');
filelist=filelist(1:size(filelist));
filenames={'<New File>'};
for i=1:size(filelist,1)
    filenames=[filenames ; filelist(i).name];
end
set(handles.Filename,'String',filenames);
contents=get(handles.Filename,'String'); %Should default to <New File>
% Set some UI stuff
set(handles.mgridfile,'Enable','off');
set(handles.table,'Data',handles.ns_array);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes VMECedit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VMECedit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Filename.
function Filename_Callback(hObject, eventdata, handles)
% hObject    handle to Filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Filename contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Filename
contents = cellstr(get(hObject,'String'));
filename=contents{get(hObject,'Value')};
if get(hObject,'Value') > 1
    handles.data=read_vmec_input(filename);
end
% Must update non-data values
handles.ns_array=handles.data.ns_array;
handles.ftol_array=handles.data.ftol_array;
handles.am=handles.data.am;
handles.ai=handles.data.ai;
handles.ac=handles.data.ac;
handles.raxis=handles.data.raxis;
handles.zaxis=handles.data.zaxis;
handles.extcurr=handles.data.extcur;
handles.rbc=handles.data.rbc;
handles.zbs=handles.data.zbs;
handles.nfp_index=handles.data.nfp;
handles.ntor_index=handles.data.ntor;
handles.mpol_index=handles.data.mpol;
% Now Update the GUI
set(handles.nfp,'String',num2str(handles.nfp_index));
set(handles.mpol,'String',num2str(handles.mpol_index));
set(handles.ntor,'String',num2str(handles.ntor_index));
set(handles.delt,'String',num2str(handles.data.delt));
set(handles.niter,'String',num2str(handles.data.niter));
set(handles.nstep,'String',num2str(handles.data.nstep));
set(handles.nvacskip,'String',num2str(handles.data.nvacskip));
set(handles.gamma,'String',num2str(handles.data.gamma));
set(handles.phiedge,'String',num2str(handles.data.phiedge));
set(handles.curtor,'String',num2str(handles.data.curtor));
set(handles.ncurrmenu,'Value',round(handles.data.ncurr+1));
if handles.data.lfreeb
    set(handles.boundary,'Value',1);
    set(handles.mgridfile,'String',handles.data.mgrid_file);
    set(handles.mgridfile,'Enable','on');
    set(handles.nzeta,'String',num2str(handles.data.nzeta));
    set(handles.nzeta,'Enable','on');
    set(handles.nextcurr,'String',num2str(numel(handles.data.extcur)));
    set(handles.nextcurr,'Enable','on');
else
    set(handles.boundary,'Value',0);
    set(handles.mgridfile,'String','mgrid_file');
    set(handles.mgridfile,'Enable','off');
    set(handles.nzeta,'String',num2str(handles.data.nzeta));
    set(handles.nzeta,'Enable','off');
    set(handles.nextcurr,'String',num2str(numel(handles.data.extcur)));
    set(handles.nextcurr,'Enable','off');
end

% Now update the handles object
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nfp_Callback(hObject, eventdata, handles)
% hObject    handle to nfp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nfp as text
%        str2double(get(hObject,'String')) returns contents of nfp as a double
handles.nfp_index=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function nfp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nfp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mpol_Callback(hObject, eventdata, handles)
% hObject    handle to mpol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mpol as text
%        str2double(get(hObject,'String')) returns contents of mpol as a double

handles.mpol_index=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function mpol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mpol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ntor_Callback(hObject, eventdata, handles)
% hObject    handle to ntor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ntor as text
%        str2double(get(hObject,'String')) returns contents of ntor as a
%        double
handles.npol_index=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function ntor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ntor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in boundary.
function boundary_Callback(hObject, eventdata, handles)
% hObject    handle to boundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.ncurrmenu,'String'));
stemp=contents{get(handles.ncurrmenu,'Value')};
selection = str2double(stemp(1:2));
boundary=get(handles.boundary,'Value');
if boundary == 1
    set(handles.mgridfile,'Enable','on');
    set(handles.nzeta,'Enable','on');
    set(handles.nextcurr,'Enable','on');
    if selection==0
        set(handles.array_menu,'String',{'NS_ARRAY' 'FTOL_ARRAY' 'AM' 'AI' 'RAXIS'...
            'ZAXIS' 'EXTCURR' 'RBC' 'ZBS'});
    elseif selection==1
        set(handles.array_menu,'String',{'NS_ARRAY' 'FTOL_ARRAY' 'AM' 'AC' 'RAXIS'...
            'ZAXIS' 'EXTCURR' 'RBC' 'ZBS'});
    end
else
    set(handles.mgridfile,'Enable','off');
    set(handles.nzeta,'Enable','off');
    set(handles.nextcurr,'Enable','off');
    if selection==0
        set(handles.array_menu,'String',{'NS_ARRAY' 'FTOL_ARRAY' 'AM' 'AI' 'RAXIS'...
            'ZAXIS' 'RBC' 'ZBS'});
    elseif selection==1
        set(handles.array_menu,'String',{'NS_ARRAY' 'FTOL_ARRAY' 'AM' 'AC' 'RAXIS'...
            'ZAXIS' 'RBC' 'ZBS'});
    end
end
    
% Hint: get(hObject,'Value') returns toggle state of boundary



function mgridfile_Callback(hObject, eventdata, handles)
% hObject    handle to mgridfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mgridfile as text
%        str2double(get(hObject,'String')) returns contents of mgridfile as a double


% --- Executes during object creation, after setting all properties.
function mgridfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mgridfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in table.
function table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
contents=cellstr(get(handles.array_menu,'String'));
index=eventdata.Indices;
val=str2double(eventdata.EditData);
switch contents{get(handles.array_menu,'Value')}
    case 'NS_ARRAY'
        temp=handles.ns_array;
        if index(1)>size(temp,1)
            temp=zeros(index(1),size(temp,2));
        end
        if index(2)>size(temp,2)
            temp=zeros(size(temp,1),index(2));
        end
        temp(1:size(handles.ns_array,1),1:size(handles.ns_array,2))=handles.ns_array;
        handles.ns_array=temp;
        handles.ns_array(index(1),index(2))=val;
    case 'FTOL_ARRAY'
        temp=handles.ftol_array;
        if index(1)>size(temp,1)
            temp=zeros(index(1),size(temp,2));
        end
        if index(2)>size(temp,2)
            temp=zeros(size(temp,1),index(2));
        end
        temp(1:size(handles.ftol_array,1),1:size(handles.ftol_array,2))=handles.ftol_array;
        handles.ftol_array=temp;
        handles.ftol_array(index(1),index(2))=val;
    case 'AM'
        if index(2) < 11
            handles.am(index(1),index(2))=val;
        end
        r=[0:.01:1];
        func=ones(1,size(r,2))*handles.am(1);
        for i=2:size(handles.am,2)
            func(:)=func(:)+handles.am(i).*r(:).^(i-1);
        end
        plot(handles.plot,r,func);
    case 'AI'
        if index(2) < 11
            handles.ai(index(1),index(2))=val;
        end
        r=[0:.01:1];
        func=ones(1,size(r,2))*handles.ai(1);
        for i=2:size(handles.ai,2)
            func(:)=func(:)+handles.ai(i).*r(:).^(i-1);
        end
        plot(handles.plot,r,func);
    case 'AC'
        if index(2) < 11
            handles.ac(index(1),index(2))=val;
        end
        r=[0:.01:1];
        func=ones(1,size(r,2))*handles.ac(1);
        for i=2:size(handles.ac,2)
            func(:)=func(:)+handles.ac(i).*r(:).^(i-1);
        end
        plot(handles.plot,r,func);
    case 'EXTCURR'
        temp=handles.extcurr;
        if index(1)>size(temp,1)
            temp=zeros(index(1),size(temp,2));
        end
        if index(2)>size(temp,2)
            temp=zeros(size(temp,1),index(2));
        end
        temp(1:size(handles.extcurr,1),1:size(handles.extcurr,2))=handles.extcurr;
        handles.extcurr=temp;
        handles.extcurr(index(1),index(2))=val;
    case 'RAXIS'
        handles.raxis(index(1),index(2))=val;
        zeta=[0:2*pi/360:pi];
        r=zeros(1,size(zeta,2));
        z2=zeros(1,size(zeta,2));
        for i=1:size(handles.raxis,2)
            r(:)=r(:)+handles.raxis(i).*cos(-(i-1).*zeta(:)*handles.nfp_index);
            z2(:)=z2(:)+handles.zaxis(i).*sin(-(i-1).*zeta(:)*handles.nfp_index);
        end
        theta=0:2*pi/360.:2*pi;
        x=zeros(size(theta,2),size(zeta,2));
        y=zeros(size(theta,2),size(zeta,2));
        z=zeros(size(theta,2),size(zeta,2));
        for i=1:size(theta,2)
                x(i,:)=r(:).*cos(zeta(:));
                y(i,:)=r(:).*sin(zeta(:));
                z(i,:)=z2(:);
        end
        plot3(x,y,z,'k')
        rotate3d on;
        axis equal
    case 'ZAXIS'
        handles.zaxis(index(1),index(2))=val;
        zeta=[0:2*pi/360:pi];
        r=zeros(1,size(zeta,2));
        z2=zeros(1,size(zeta,2));
        for i=1:size(handles.raxis,2)
            r(:)=r(:)+handles.raxis(i).*cos(-(i-1).*zeta(:)*handles.nfp_index);
            z2(:)=z2(:)+handles.zaxis(i).*sin(-(i-1).*zeta(:)*handles.nfp_index);
        end
        theta=0:2*pi/360.:2*pi;
        x=zeros(size(theta,2),size(zeta,2));
        y=zeros(size(theta,2),size(zeta,2));
        z=zeros(size(theta,2),size(zeta,2));
        for i=1:size(theta,2)
                x(i,:)=r(:).*cos(zeta(:));
                y(i,:)=r(:).*sin(zeta(:));
                z(i,:)=z2(:);
        end
        plot3(x,y,z,'k')
        rotate3d on;
        axis equal
    case 'RBC'
        handles.rbc(index(1),index(2))=val;
        n=35;
        m=71;
        theta=0:2*pi/m:2*pi;
        zeta=0:pi/n/handles.nfp_index:pi/handles.nfp_index;
        r=cfunct_old(theta,zeta,handles.rbc',handles.nfp_index);
        z=sfunct_old(theta,zeta,handles.zbs',handles.nfp_index);
        lambda=zeros([1 m+1 n+1]);
        cla;
        isotoro(r,z,zeta,1);
        camlight headlight;
        rotate3d on;
    case 'ZBS'
        handles.zbs(index(1),index(2))=val;
        n=35;
        m=71;
        theta=0:2*pi/m:2*pi;
        zeta=0:pi/n/handles.nfp_index:pi/handles.nfp_index;
        r=cfunct_old(theta,zeta,handles.rbc',handles.nfp_index);
        z=sfunct_old(theta,zeta,handles.zbs',handles.nfp_index);
        lambda=zeros([1 m+1 n+1]);
        cla;
        isotoro(r,z,zeta,1);
        camlight headlight;
        rotate3d on;
end
guidata(hObject, handles);


% --- Executes on selection change in array_menu.
function array_menu_Callback(hObject, eventdata, handles)
% hObject    handle to array_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns array_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from array_menu
contents=cellstr(get(hObject,'String'));

%set(handles.rowlabel,'Visible','off');
set(handles.table,'RowName',[1]);
set(handles.collabel,'String','Index');
set(handles.table,'ColumnName',0:10);
rotate3d off;
switch contents{get(hObject,'Value')}
    case 'NS_ARRAY'
        set(handles.table,'Data',handles.ns_array);
        set(handles.table,'ColumnEditable',true(1,size(handles.ns_array,2)));
    case 'FTOL_ARRAY'
        set(handles.table,'Data',handles.ftol_array);
        set(handles.table,'ColumnEditable',true(1,size(handles.ftol_array,2)));
    case 'AM'
        set(handles.table,'Data',handles.am);
        set(handles.table,'ColumnEditable',true(1,size(handles.am,2)));
        r=[0:.01:1];
        func=ones(1,size(r,2))*handles.am(1);
        for i=2:size(handles.am,2)
            func(:)=func(:)+handles.am(i).*r(:).^(i-1);
        end
        plot(handles.plot,r,func);
    case 'AI'
        set(handles.table,'Data',handles.ai);
        set(handles.table,'ColumnEditable',true(1,size(handles.ai,2)));
        r=[0:.01:1];
        func=ones(1,size(r,2))*handles.ai(1);
        for i=2:size(handles.ai,2)
            func(:)=func(:)+handles.ai(i).*r(:).^(i-1);
        end
        plot(handles.plot,r,func);
    case 'AC'
        set(handles.table,'Data',handles.ac);
        set(handles.table,'ColumnEditable',true(1,size(handles.ac,2)));
        r=[0:.01:1];
        func=ones(1,size(r,2))*handles.ac(1);
        for i=2:size(handles.ac,2)
            func(:)=func(:)+handles.ac(i).*r(:).^(i-1);
        end
        plot(handles.plot,r,func);
    case 'EXTCURR'
        set(handles.table,'Data',handles.extcurr);
        set(handles.table,'ColumnEditable',true(1,size(handles.extcurr,2)));
    case 'RAXIS'
        set(handles.table,'Data',handles.raxis);
        set(handles.table,'ColumnEditable',true(1,size(handles.raxis,2)));
        zeta=[0:2*pi/360:pi];
        r=zeros(1,size(zeta,2));
        z2=zeros(1,size(zeta,2));
        for i=1:size(handles.raxis,2)
            r(:)=r(:)+handles.raxis(i).*cos(-(i-1).*zeta(:)*handles.nfp_index);
            z2(:)=z2(:)+handles.zaxis(i).*sin(-(i-1).*zeta(:)*handles.nfp_index);
        end
        theta=0:2*pi/360.:2*pi;
        x=zeros(size(theta,2),size(zeta,2));
        y=zeros(size(theta,2),size(zeta,2));
        z=zeros(size(theta,2),size(zeta,2));
        for i=1:size(theta,2)
                x(i,:)=r(:).*cos(zeta(:));
                y(i,:)=r(:).*sin(zeta(:));
                z(i,:)=z2(:);
        end
        plot3(x,y,z,'k')
        rotate3d on;
        axis equal
    case 'ZAXIS'
        set(handles.table,'Data',handles.zaxis);
        set(handles.table,'ColumnEditable',true(1,size(handles.zaxis,2)));
        zeta=[0:2*pi/360:pi];
        r=zeros(1,size(zeta,2));
        z2=zeros(1,size(zeta,2));
        for i=1:size(handles.raxis,2)
            r(:)=r(:)+handles.raxis(i).*cos(-(i-1).*zeta(:)*handles.nfp_index);
            z2(:)=z2(:)+handles.zaxis(i).*sin(-(i-1).*zeta(:)*handles.nfp_index);
        end
        theta=0:2*pi/360.:2*pi;
        x=zeros(size(theta,2),size(zeta,2));
        y=zeros(size(theta,2),size(zeta,2));
        z=zeros(size(theta,2),size(zeta,2));
        for i=1:size(theta,2)
                x(i,:)=r(:).*cos(zeta(:));
                y(i,:)=r(:).*sin(zeta(:));
                z(i,:)=z2(:);
        end
        plot3(x,y,z,'k')
        rotate3d on;
        axis equal
    case 'RBC'
        handles.ntor_index=str2double(get(handles.ntor,'String'));
        rownames=-handles.ntor_index:1:handles.ntor_index;
        handles.mpol_index=str2double(get(handles.mpol,'String'));
        colnames=0:1:handles.mpol_index;
        set(handles.table,'RowName',rownames);
        set(handles.table,'ColumnName',colnames);
%        set(handles.rowlabel,'Visible','on');
        set(handles.collabel,'String','mpol');
        set(handles.table,'Data',handles.rbc);
        set(handles.table,'ColumnEditable',true(1,size(handles.rbc,2)));
        n=35;
        m=71;
        theta=0:2*pi/m:2*pi;
        zeta=0:pi/n/handles.nfp_index:pi/handles.nfp_index;
        r=cfunct_old(theta,zeta,handles.rbc',handles.nfp_index);
        z=sfunct_old(theta,zeta,handles.zbs',handles.nfp_index);
        lambda=zeros([1 m+1 n+1]);
        cla;
        isotoro(r,z,zeta,1);
        camlight headlight;
        rotate3d on;
    case 'ZBS'
        handles.ntor_index=str2double(get(handles.ntor,'String'));
        rownames=-handles.ntor_index:1:handles.ntor_index;
        handles.mpol_index=str2double(get(handles.mpol,'String'));
        colnames=0:1:handles.mpol_index;
        set(handles.table,'RowName',rownames);
        set(handles.table,'ColumnName',colnames);
        %set(handles.rowlabel,'Visible','on');
        set(handles.collabel,'String','mpol');
        set(handles.table,'Data',handles.zbs);
        set(handles.table,'ColumnEditable',true(1,size(handles.zbs,2)));
        n=35;
        m=71;
        theta=0:2*pi/m:2*pi;
        zeta=0:pi/n/handles.nfp_index:pi/handles.nfp_index;
        r=cfunct_old(theta,zeta,handles.rbc',handles.nfp_index);
        z=sfunct_old(theta,zeta,handles.zbs',handles.nfp_index);
        lambda=zeros([1 m+1 n+1]);
        cla;
        isotoro(r,z,zeta,1);
        camlight headlight;
        rotate3d on;
end

% --- Executes during object creation, after setting all properties.
function array_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to array_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ncurrmenu.
function ncurrmenu_Callback(hObject, eventdata, handles)
% hObject    handle to ncurrmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ncurrmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ncurrmenu
contents = cellstr(get(handles.ncurrmenu,'String'));
stemp=contents{get(handles.ncurrmenu,'Value')};
selection = str2double(stemp(1:2));
boundary=get(handles.boundary,'Value');
if boundary == 1
    set(handles.mgridfile,'Enable','on');
    set(handles.nzeta,'Enable','on');
    if selection==0
        set(handles.array_menu,'String',{'NS_ARRAY' 'FTOL_ARRAY' 'AM' 'AI' 'RAXIS'...
            'ZAXIS' 'EXTCURR' 'RBC' 'ZBS'});
    elseif selection==1
        set(handles.array_menu,'String',{'NS_ARRAY' 'FTOL_ARRAY' 'AM' 'AC' 'RAXIS'...
            'ZAXIS' 'EXTCURR' 'RBC' 'ZBS'});
    end
else
    set(handles.mgridfile,'Enable','off');
    set(handles.nzeta,'Enable','off');
    if selection==0
        set(handles.array_menu,'String',{'NS_ARRAY' 'FTOL_ARRAY' 'AM' 'AI' 'RAXIS'...
            'ZAXIS' 'RBC' 'ZBS'});
    elseif selection==1
        set(handles.array_menu,'String',{'NS_ARRAY' 'FTOL_ARRAY' 'AM' 'AC' 'RAXIS'...
            'ZAXIS' 'RBC' 'ZBS'});
    end
end


% --- Executes during object creation, after setting all properties.
function ncurrmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ncurrmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delt_Callback(hObject, eventdata, handles)
% hObject    handle to delt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delt as text
%        str2double(get(hObject,'String')) returns contents of delt as a double


% --- Executes during object creation, after setting all properties.
function delt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nzeta_Callback(hObject, eventdata, handles)
% hObject    handle to nzeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nzeta as text
%        str2double(get(hObject,'String')) returns contents of nzeta as a double


% --- Executes during object creation, after setting all properties.
function nzeta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nzeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function niter_Callback(hObject, eventdata, handles)
% hObject    handle to niter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of niter as text
%        str2double(get(hObject,'String')) returns contents of niter as a double


% --- Executes during object creation, after setting all properties.
function niter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to niter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nstep_Callback(hObject, eventdata, handles)
% hObject    handle to nstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nstep as text
%        str2double(get(hObject,'String')) returns contents of nstep as a double


% --- Executes during object creation, after setting all properties.
function nstep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nvacskip_Callback(hObject, eventdata, handles)
% hObject    handle to nvacskip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nvacskip as text
%        str2double(get(hObject,'String')) returns contents of nvacskip as a double


% --- Executes during object creation, after setting all properties.
function nvacskip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nvacskip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gamma_Callback(hObject, eventdata, handles)
% hObject    handle to gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gamma as text
%        str2double(get(hObject,'String')) returns contents of gamma as a double


% --- Executes during object creation, after setting all properties.
function gamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function phiedge_Callback(hObject, eventdata, handles)
% hObject    handle to phiedge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phiedge as text
%        str2double(get(hObject,'String')) returns contents of phiedge as a double


% --- Executes during object creation, after setting all properties.
function phiedge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phiedge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function curtor_Callback(hObject, eventdata, handles)
% hObject    handle to curtor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of curtor as text
%        str2double(get(hObject,'String')) returns contents of curtor as a double


% --- Executes during object creation, after setting all properties.
function curtor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to curtor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outfilename_Callback(hObject, eventdata, handles)
% hObject    handle to outfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outfilename as text
%        str2double(get(hObject,'String')) returns contents of outfilename as a double


% --- Executes during object creation, after setting all properties.
function outfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in makefile.
function makefile_Callback(hObject, eventdata, handles)
% hObject    handle to makefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=get(handles.outfilename,'String');
fid=fopen(filename,'wt');
fprintf(fid,'%s\n',['!----- Created by VMECedit ' datestr(now) ' -----']);
fprintf(fid,'%s\n','&INDATA');
fprintf(fid,'%s\n','!----- Runtime Parameters -----');
write_namelist_flt(fid,'DELT',str2double(get(handles.delt,'String')));
write_namelist_int(fid,'NITER',str2double(get(handles.niter,'String')));
write_namelist_int(fid,'NSTEP',str2double(get(handles.nstep,'String')));
write_namelist_flt(fid,'TCON0',1.0);
write_namelist_vec(fid,'NS_ARRAY',handles.ns_array,'int');
write_namelist_vec(fid,'FTOL_ARRAY',handles.ftol_array);
fprintf(fid,'%s\n','!----- Grid Parameters -----');
fprintf(fid,'%s\n','  LASYM = F');
write_namelist_int(fid,'NFP',str2double(get(handles.nfp,'String')));
write_namelist_int(fid,'MPOL',str2double(get(handles.mpol,'String')));
write_namelist_int(fid,'NTOR',str2double(get(handles.ntor,'String')));
write_namelist_flt(fid,'PHIEDGE',str2double(get(handles.phiedge,'String')));
fprintf(fid,'%s\n','!----- Free Boundary Parameters -----');
if get(handles.boundary,'Value')
    fprintf(fid,'%s\n','  LFREEB = T');
    write_namelist_str(fid,'MGRID_FILE',get(handles.mgridfile,'String'));
    write_namelist_int(fid,'NTHETA',str2double(get(handles.ntor,'String'))*2+6);
    write_namelist_int(fid,'NZETA',str2double(get(handles.nzeta,'String')));
    write_namelist_vec(fid,'EXTCUR',handles.extcurr,'int');
    write_namelist_int(fid,'NVACSKIP',str2double(get(handles.nvacskip,'String')));
else
    fprintf(fid,'%s\n','  LFREEB = F');
    fprintf(fid,'%s = ''%s''\n','  MGRID_FILE','NONE');
end
fprintf(fid,'%s\n','!----- Pressure Parameters -----');
write_namelist_flt(fid,'GAMMA',str2double(get(handles.gamma,'String')));
write_namelist_flt(fid,'BLOAT',1.0);
write_namelist_flt(fid,'SPRES_PED',1.0);
write_namelist_vec(fid,'AM',handles.am);
fprintf(fid,'%s\n','!----- Current/Iota Parameters -----');
write_namelist_flt(fid,'CURTOR',str2double(get(handles.curtor,'String')));
write_namelist_int(fid,'NCURR',get(handles.ncurrmenu,'Value')-1);
write_namelist_int(fid,'AC_FORM',0);
write_namelist_vec(fid,'AI',handles.ai);
write_namelist_vec(fid,'AC',handles.ac);
fprintf(fid,'%s\n','!----- Axis Parameters -----');
write_namelist_vec(fid,'RAXIS',handles.raxis);
write_namelist_vec(fid,'ZAXIS',handles.zaxis);
fprintf(fid,'%s\n','!----- Boundary Parameters -----');
write_namelist_arr(fid,'RBC',handles.rbc,...
    -handles.ntor_index,handles.ntor_index,0,handles.mpol_index);
write_namelist_arr(fid,'ZBS',handles.zbs,...
    -handles.ntor_index,handles.ntor_index,0,handles.mpol_index);
fprintf(fid,'%s\n\n','/');
fclose(fid);



function nextcurr_Callback(hObject, eventdata, handles)
% hObject    handle to nextcurr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nextcurr as text
%        str2double(get(hObject,'String')) returns contents of nextcurr as a double
temp=handles.extcurr;
handles.extcurr=zeros(1,str2double(get(hObject,'String')));
if size(temp,2) >= str2double(get(hObject,'String'))
    handles.extcurr(:)=temp(1:str2double(get(hObject,'String')));
else
    handles.extcurr(1:size(temp,2))=temp(:);
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function nextcurr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nextcurr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
