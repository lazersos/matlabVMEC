function varargout = VMECplot(varargin)
%VMECPLOT([data]) GUI for interactive VMEC plots.
%   VMECPLOT is a MATLAB Graphical User Interface (GUI) for
%   visualizing the data stored in the VMEC output file (wout.*).  It
%   allows the user to plot various quantities in the output file in
%   various ways.  Command line data read may be passed to VMECplot for
%   plotting.
%
%   Options:
%
%       VMECplot(vmec_data) Plots data in VMEC_DATA as read by read_vmec.
%
%       VMECplot(ves_data) Vessel data added to plots as read by
%       read_vessel.
%
%       VMECplot('nu',nu) Defalut the number of poloidal gridpoints to nu.
%       This overrides an internal calculation of nu based on mpol.
%
%       VMECplot('nv',nv) Default the number of toroidal gridpoints to nv.
%       This overrides an internal calculation of nv based on ntor.  This
%       is ignored if ntor=0 (Tokamak).
%
%   See also read_vmec, read_vessel, sfunct, cfunct, torocont, isotoro,
%   plot_vessel.
%
%   Example:
%      data=read_vmec('wout.test');
%      VMECplot(data)
%      ves_data=read_vessel('vessel.dat');
%      VMECplot(data,ves_data)
%      VMECplot(ves_data)
%
%   Written by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:    2.5
%   Date:       03/17/2012

% SAl-10/20/10
% Accepts data as a command line argument for plotting.
%
% SAL-10/14/10
% Edited to Version 1.0 
%
% SAL-10/7/10
% Poloidal 2D plots now plot flux surface vs. zeta
% Toroidal slider now only works over a field period
% Poloidal 2D plots now only plot over a field period in zeta
% Flux Surface 2D plots now properly added surfaces via slider
% Fluxcut 2D plots not only plot over a field period in zeta
% Saveas function now outputs to all available extensions in saveas
%      ai,bmp,emf,eps,fig,jpg,m,pbm,pcx,pdf,pgm,png,ppm,tif
%
% SAL-12/22/10
% Fixed an issue when trying to plot files with ntor=0 (tokamaks) thanks to
% Peter Kavran for bringing this to my attention.
%
% SAL-01/10/11
% Will now accept vaccum vessel data for plotting
%    VMECplot(vac_data)
% Also fixed issue with 3D plotting of surfaces where -zeta was used.  Now
% uses zeta
%
% SAL-01/31/11
% All values now assumed to be on the full mesh. 1D values plotted in flux
% space.
%
% SAL-02/01/11
% Updated 1D values.
% Dynamic menus.
% Uses VMEC nu represenation (lambda)
% Flux Sufaces and LPK plots now have magnetic axes
% Removed Field Line plots as it was misleading
%
% SAL-03/11/11
% Numerous ploting bug fixes.
% Added command line option for PIES plotting support, override of nu and
% nv calculations.
% Also fixed some bugs plotting Tokamaks.
%
% SAL-03/31/11
% Added 'diag' option to allow diagnostic plots (force norm and fourier)
%
% SAL-04/06/11
% Added ability to plot boozer data from booz_xform
%
% SAL-10/17/11
% Switched to faster vectorized form of the fourier transform function.
%
% SAL-02/23/12
% Added support for SPEC code
%
% SAL-03/17/12
% Fixed calculation of J (needed to be divided by g)
% Fixed labeling of curru and currv plots
% Fixed non-stellarator symmetric calculations
% Fixed error when attempting to isotoro plot the axis (now plots 3D axis)
% Updated version to 2.1
%
% SAL-03/22/12
% Fixed file output with multiple periods in filename.
% LPK plot now available for nfp=1 ntor>0 equilibria.
%
% SAL-04/19/12
% Added support for SPEC convergence and grid plot
%
% SAL-08/27/13
% Added support for NSTAB equilibria
%
% SAL-12/21/13
% Added support for VMEC Flow variable plotting
%
% SAL-03/07/20
% Added support for SIESTA Plotting

% Edit the above text to modify the response to help VMECplot

% Last Modified by GUIDE v2.5 06-Aug-2010 15:40:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @VMECplot_OpeningFcn, ...
                   'gui_OutputFcn',  @VMECplot_OutputFcn, ...
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


% --- Executes just before VMECplot is made visible.
function VMECplot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VMECplot (see VARARGIN)

numdefargs=3;
% Choose default command line output for VMECplot
handles.output = hObject;
% Get the filelist
filelist=dir('wout*');
filelist=filelist(1:size(filelist));
% Set some initial stuff
handles.cuttype='text';
handles.rval=1;
handles.thetaval=1;
handles.zetaval=1;
handles.plotves=0;
handles.nuoverride=0;
handles.nvoverride=0;
handles.diag_plt=0;
% Handle no wout files
if isempty(filelist) && (nargin ==3)
    disp(' - No wout files found in current directory.');
    axis([0 1 0 1]);
    cla;
    text(0.05,0.6,'!!!!No wout files found!!!!!',...
        'Color','red','FontSize',30);
    set(handles.statustext,'String',...
        'NO WOUT FILES');
    pause(0.01);
    zoom off
    colorbar('off')
    rotate3d off
    set(handles.rslide,'Enable','off');
    set(handles.theslide,'Enable','off');
    set(handles.torslide,'Enable','off');
    set(handles.rcut,'Enable','off');
    set(handles.thetacut,'Enable','off');
    set(handles.zetacut,'Enable','off');
    set(handles.fluxcut,'Enable','off');
    set(handles.polcut,'Enable','off');
    set(handles.rzcut,'Enable','off');
    set(handles.threedcut,'Enable','off');
    set(handles.filename,'Enable','off');
    set(handles.plottype,'Enable','off');
    guiupdate(handles);
    return
end
guidata(hObject, handles);
if (nargin == numdefargs)
    % Set the Filenames
    set(handles.filename,'String',{filelist.name});
    contents=get(handles.filename,'String');
    % Read the currently chosen file
    handles.data=read_vmec(contents{get(handles.filename,'Value')});
elseif (nargin > numdefargs)
    i=1;
    while i<=nargin-numdefargs
        if isstruct(varargin{i})
            switch varargin{i}.datatype
                case {'wout','pies_out','boozer','SPEC','FIELDLINES','nout','siesta','BEAMS3D'}
                    handles.data=varargin{i};
                    if isempty(filelist)
                        set(handles.filename,'String',{'User Data'});
                    else
                        set(handles.filename,'String',{'User Data' filelist.name});
                    end
                case 'vessel'
                    handles.plotves=1;
                    handles.ves_data=varargin{i};
                    disp(' - Vessel Data Detected');
                otherwise
                    disp(['ERROR: Unknown datatype:' varargin{i}.datatype]);
            end
        else
            switch varargin{i}
                case 'nu'
                    i=i+1;
                    handles.mpol=varargin{i};
                    handles.nuoverride=1;
                case 'nv'
                    i=i+1;
                    handles.ntor=varargin{i};
                    handles.nvoverride=1;
                case 'diag'
                    handles.diag_plt=1;
                otherwise
                    disp(['Error: Unknown option:' varargin{i}]);
            end
        end
        i=i+1;
    end
    % Read in a datafile if none was supplied
    if ~isfield(handles,'data')
        set(handles.filename,'String',{filelist.name});
        contents=get(handles.filename,'String');
        % Read the currently chosen file
        handles.data=read_vmec(contents{get(handles.filename,'Value')});
    end
end
% Handle a VMEC Runtime Error File
if isfield(handles.data,'ierr_vmec')
    if (handles.data.ierr_vmec && (handles.data.ierr_vmec ~= 4))
        disp(strcat(' - VMECplot has detected an error in :',...
            contents{get(handles.filename,'Value')}));
        set(handles.statustext,'String',...
            'VMEC ERROR');
        axis([0 1 0 1]);
        cla;
        text(0.05,0.6,'VMEC Runtime Error Detected',...
            'Color','red','FontSize',24);
        text(0.25,0.5,...
            strcat('File: ',contents{get(handles.filename,'Value')}),...
            'Interpreter','none');
        text(0.25,0.4,strcat('ierr_vmec=',num2str(handles.data.ierr_vmec)),...
            'Interpreter','none');
        pause(.01);
        zoom off
        colorbar('off')
        rotate3d off
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','off');
        set(handles.rcut,'Enable','off');
        set(handles.thetacut,'Enable','off');
        set(handles.zetacut,'Enable','off');
        set(handles.fluxcut,'Enable','off');
        set(handles.polcut,'Enable','off');
        set(handles.rzcut,'Enable','off');
        set(handles.threedcut,'Enable','off');
        return
    end
end
% Quick handle of BEAMS3D
if ~isfield(handles.data,'ntor')
    handles.data.ntor = handles.data.nphi;
end
if ~isfield(handles.data,'nfp')
    handles.data.nfp = 1;
end
if ~isfield(handles.data,'ns')
    handles.data.ns = 128;
end
% Set up grids
nsin=5.;
if handles.data.ntor == 0
    handles.ntor=32;
    if ~handles.nuoverride, handles.mpol=64; end
    handles.theta=0:2*pi/(handles.mpol-1):2*pi;
    handles.zeta=0:2*pi/(handles.ntor-1):2*pi;
    disp(strcat(' - ns=',num2str(handles.data.ns)));
    disp(strcat(' - ntheta=',num2str(handles.mpol)));
    disp(strcat(' - nzeta=',num2str(handles.ntor)));
    disp(' -- Loading Data, Please wait --');
    set(handles.torslide,'Enable','off');
    set(handles.torslide,'Max',1.0);
    set(handles.torslide,'Min',1.0);
    set(handles.torslide,'SliderStep',[1.0 1.0]);
    set(handles.torslide,'Value',1.0);
else
    if ~handles.nvoverride
        handles.ntor=max(4*handles.data.ntor+1,32);
        %handles.ntor=max(nsin*(handles.data.ntor+1)*handles.data.nfp,32);
    end
    handles.mpol=64;
    %if ~handles.nuoverride, handles.mpol=handles.data.nu;end
    disp(strcat(' - ns=',num2str(handles.data.ns)));
    disp(strcat(' - ntheta=',num2str(handles.mpol)));
    disp(strcat(' - nzeta=',num2str(handles.ntor)));
    disp(' -- Loading Data, Please wait --');
    handles.zeta=0:2*pi/double(handles.ntor-1):2*pi;          %nzeta and ntor+1 elements
    handles.theta=0:2*pi/double(handles.mpol-1):2*pi;
    set(handles.torslide,'Max',handles.ntor);
    set(handles.torslide,'Min',1.0);
    set(handles.torslide,'SliderStep',[1.0 1.0]./double(int32(handles.ntor)-1));
    set(handles.torslide,'Value',1.0);
    set(handles.torslide,'Enable','off');
end
% Set Defaults
set(handles.torslide,'Max',handles.ntor);
set(handles.torslide,'Min',1.0);
set(handles.torslide,'SliderStep',[1.0 1.0]./double(int32(handles.ntor)-1));
set(handles.torslide,'Value',1.0);
%set(handles.torslide,'Max',round(handles.ntor/handles.data.nfp)+1);
%set(handles.torslide,'Min',1.0);
%set(handles.torslide,'SliderStep',[1.0 1.0]./double(round(int32(handles.ntor)/int32(handles.data.nfp))+1-1));
%set(handles.torslide,'Value',1.0);
set(handles.theslide,'Max',handles.mpol);
set(handles.theslide,'Min',1.0);
set(handles.theslide,'SliderStep',[1.0 1.0]./double(int32(handles.mpol)-1));
set(handles.theslide,'Value',1.0);
set(handles.rslide,'Max',handles.data.ns);
set(handles.rslide,'Min',1.0);
set(handles.rslide,'SliderStep',[1.0 1.0]./double(int32(handles.data.ns)-1));
set(handles.rslide,'Value',1.0);
set(handles.rslide,'Enable','off');
set(handles.theslide,'Enable','off');
set(handles.rtext,'String','Flux Surf');
set(handles.thetext,'String','Theta');
set(handles.tortext,'String','Zeta');
set(handles.rcut,'Enable','off');
set(handles.rcut,'Value',1.0);
set(handles.thetacut,'Enable','off');
set(handles.zetacut,'Enable','off');
set(handles.fluxcut,'Enable','off');
set(handles.polcut,'Enable','off');
set(handles.rzcut,'Enable','off');
set(handles.threedcut,'Enable','off');
set(handles.cutplane,'SelectedObject',handles.rcut);
% Update handles structure
guidata(hObject, handles);
% Transform Spectrum
handles=transf(hObject,handles);
% Setup the Available plot types
vfields={'iotaf' 'presf' 'Dmerc' ...
    'Dshear' 'Dwell' 'Dcurr' 'Dgeod' 'jdotb' 'bdotgradv' 'beta_vol' ...
    'phip' 'buco' 'bvco' 'phi' 'vp' 'overr' 'jcuru' 'jcurv' 'specw' ...
    'dpdr' 'resid' 'nisla' 'correc' 'islw' ...
    'correc2' 'jdev' 'edgeo' 'iota_step' 'fits','errp','mu',...
    'pmap','omega','tpotb'};
mfields={'b' 'lam' 'p' 'g' 'b_s' 'b_u' 'b_v' 'bs' 'bu' 'bv' 'br' 'bphi' 'bz'...
    'currs' 'curru' 'currv' 'jr' 'jphi' 'jz' 'brho' 'bnorm' 'btheta' 'jpara' 'press'...
    'fnx' 'fny' 'fnz' 'fn' 'rbc' 'rbs' 'zbc' 'zbs' 'jxbx' 'jxby' 'jxbz' ...
    'dpds' 'ptran' 'prot' 'prpr' 'U_rot' 'vphi'};
optfields={};
if isfield(handles.data,'rmnc')
    optfields=[optfields 'Flux Surface' 'Flux Web'];
end
if isfield(handles,'br') && isfield(handles,'bphi') && isfield(handles,'bz')
    optfields=[optfields 'B-Field'];
end
if isfield(handles,'br_pol') && isfield(handles,'bphi') && isfield(handles,'bz_pol')
    optfields=[optfields 'B-Field(pol)'];
end
if isfield(handles,'jr') && isfield(handles,'jphi') && isfield(handles,'jz')
    optfields=[optfields 'J-Field'];
end
if handles.data.ntor > 0
    optfields=[optfields 'LPK Plot'];
end
if isfield(handles.data,'bc')
    optfields=[optfields '|B| Modes (axi)' '|B| Modes (nonaxi)'];
end
if strcmp(handles.data.datatype,'boozer')
    optfields=[optfields '|B| Field Line'];
end
if isfield(handles.data,'R_grid') && isfield(handles.data,'Z_grid')
    optfields=[optfields 'SPEC Grid'];
end
if isfield(handles.data,'R_lines') && isfield(handles.data,'Z_lines')
    optfields=[optfields 'Poincare (real)'];
end
if isfield(handles.data,'B_R')
    optfields=[optfields 'B_R'];
end
if isfield(handles.data,'B_PHI')
    optfields=[optfields 'B_PHI'];
end
if isfield(handles.data,'B_Z')
    optfields=[optfields 'B_Z'];
end
if isfield(handles.data,'S_ARR')
    optfields=[optfields 'S_ARR'];
end
if isfield(handles.data,'rho_lines') && isfield(handles.data,'th_lines')
    optfields=[optfields 'Poincare (polar)'];
    optfields=[optfields 'Poincare (rho)'];
end
temp={'Summary'; '-----1 D-----'};
for i=1:numel(vfields)
    if isfield(handles.data,vfields{i})
        temp=[temp; vfields{i}];
    end
end
temp=[temp; '-----3 D-----'];
for i=1:numel(mfields)
    if isfield(handles,mfields{i})
        temp=[temp; mfields{i}];
    end
end
temp=[temp; '-----Special-----'];
for i=1:numel(optfields)
    temp=[temp; optfields{i}];
end
set(handles.plottype,'String',temp);
set(handles.plottype,'Value',1.0);
% Update handles structure
guidata(hObject, handles);
% Update the GUI
guiupdate(handles);

% UIWAIT makes VMECplot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VMECplot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;


% --- Executes on selection change in filename.
function filename_Callback(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filename contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        filename

% Initilize some stuff
handles.cuttype='text';
handles.rval=1;
handles.thetaval=1;
handles.zetaval=1;
zoom off
colorbar('off')
rotate3d off
set(handles.rslide,'Enable','off');
set(handles.theslide,'Enable','off');
set(handles.torslide,'Enable','off');
set(handles.rcut,'Enable','off');
set(handles.thetacut,'Enable','off');
set(handles.zetacut,'Enable','off');
set(handles.fluxcut,'Enable','off');
set(handles.polcut,'Enable','off');
set(handles.rzcut,'Enable','off');
set(handles.threedcut,'Enable','off');
set(handles.cutplane,'SelectedObject',handles.rcut);
set(handles.plottype,'Value',1);
% Get the files
contents=get(hObject,'String');
handles.data=read_vmec(contents{get(hObject,'Value')});
% Handle a VMEC runtime Error
if (handles.data.ierr_vmec && (handles.data.ierr_vmec ~= 4))
    disp(strcat(' - VMECplot has detected an error in :',...
        contents{get(handles.filename,'Value')}));
    set(handles.statustext,'String',...
        'VMEC ERROR');
    axis([0 1 0 1]);
    cla;
    text(0.05,0.6,'VMEC Runtime Error Detected',...
        'Color','red','FontSize',24);
    text(0.25,0.5,...
        strcat('File: ',contents{get(handles.filename,'Value')}),...
        'Interpreter','none');
    text(0.25,0.4,strcat('ierr_vmec=',num2str(handles.data.ierr_vmec)),...
        'Interpreter','none');
    pause(.01);
    zoom off
    colorbar('off')
    rotate3d off
    set(handles.rslide,'Enable','off');
    set(handles.theslide,'Enable','off');
    set(handles.torslide,'Enable','off');
    set(handles.rcut,'Enable','off');
    set(handles.thetacut,'Enable','off');
    set(handles.zetacut,'Enable','off');
    set(handles.fluxcut,'Enable','off');
    set(handles.polcut,'Enable','off');
    set(handles.rzcut,'Enable','off');
    set(handles.threedcut,'Enable','off');
    return
end
% Setup Theta and Zeta Axes
nsin=5.;
if handles.data.ntor == 0
    handles.ntor=1;
    if ~handles.nuoverride, handles.mpol=90; end
    handles.theta=0:2*pi/(handles.mpol-1):2*pi;
    handles.zeta=0;
    disp(strcat(' - ns=',num2str(handles.data.ns)));
    disp(strcat(' - ntheta=',num2str(handles.mpol)));
    disp(strcat(' - nzeta=',num2str(handles.ntor)));
    disp(' -- Loading Data, Please wait --');
    set(handles.torslide,'Enable','off');
    set(handles.torslide,'Max',2.0);
    set(handles.torslide,'Min',1.0);
    set(handles.torslide,'SliderStep',[1.0 1.0]);
    set(handles.torslide,'Value',1.0);
else
    if ~handles.nvoverride
        %handles.ntor=nsin*(handles.data.ntor+1)*handles.data.nfp;
        handles.ntor=max(4*handles.data.ntor+1,32);
    end
    if ~handles.nuoverride, handles.mpol=90; end
    %if ~handles.nuoverride, handles.mpol=handles.data.nu; end
    disp(strcat(' - ns=',num2str(handles.data.ns)));
    disp(strcat(' - ntheta=',num2str(handles.mpol)));
    disp(strcat(' - nzeta=',num2str(handles.ntor)));
    disp(' -- Loading Data, Please wait --');
    handles.zeta=0:2*pi/double(handles.ntor-1):2*pi;          %nzeta and ntor+1 elements
    handles.theta=0:2*pi/double(handles.mpol-1):2*pi;
    set(handles.torslide,'Max',handles.ntor);
    set(handles.torslide,'Min',1.0);
    set(handles.torslide,'SliderStep',[1.0 1.0]./double(int32(handles.ntor)-1));
    set(handles.torslide,'Value',1.0);
    set(handles.torslide,'Enable','off');
end
% Setup The Sliders
set(handles.theslide,'Max',handles.mpol);
set(handles.theslide,'Min',1.0);
set(handles.theslide,'SliderStep',[1.0 1.0]./(handles.mpol-1));
set(handles.theslide,'Value',1.0);
set(handles.rslide,'Max',handles.data.ns);
set(handles.rslide,'Min',1.0);
set(handles.rslide,'SliderStep',[1.0 1.0]./(handles.data.ns-1));
set(handles.rslide,'Value',1.0);
% Now update the guidata
guidata(hObject,handles);
% Transform from Fourier to real space
handles=transf(hObject,handles);
% Setup the Available plot types
vfields={'iotaf' 'presf' 'Dmerc' ...
    'Dshear' 'Dwell' 'Dcurr' 'Dgeod' 'jdotb' 'bdotgradv' 'beta_vol' ...
    'phip' 'buco' 'bvco' 'phi' 'vp' 'overr' 'jcuru' 'jcurv' 'specw' ...
    'dpdr' 'resid' 'nisla' 'correc' 'islw' ...
    'correc2' 'jdev' 'edgeo' 'iota_step' 'fits' ...
    'pmap','omega','tpotb'};
mfields={'b' 'p' 'g' 'b_s' 'b_u' 'bs' 'b_v' 'bu' 'bv' 'curru' 'currv' 'br'...
    'bphi' 'bz' 'jr' 'jphi' 'jz' 'brho' 'bnorm' 'btheta' 'jpara' 'press'...
    'fnx' 'fny' 'fnz' 'fn' 'rbc' 'rbs' 'zbc' 'zbs' 'jxbx' 'jxby' 'jxbz' ...
    'dpds' 'ptran' 'prot' 'prpr' 'U_rot' 'vphi'};
optfields={'Flux Surface' 'Flux Web'};
if handles.data.ntor > 0
    optfields=[optfields 'LPK Plot'];
end
if isfield(handles,'br') && isfield(handles,'bphi') && isfield(handles,'bz')
    optfields=[optfields 'B-Field'];
end
if isfield(handles,'br_pol') && isfield(handles,'bphi') && isfield(handles,'bz_pol')
    optfields=[optfields 'B-Field(pol)'];
end
if isfield(handles,'jr') && isfield(handles,'jphi') && isfield(handles,'jz')
    optfields=[optfields 'J-Field'];
end
temp={'Summary'; '-----1 D-----'};
for i=1:numel(vfields)
    if isfield(handles.data,vfields{i})
        temp=[temp; vfields{i}];
    end
end
temp=[temp; '-----3 D-----'];
for i=1:numel(mfields)
    if isfield(handles,mfields{i})
        temp=[temp; mfields{i}];
    end
end
temp=[temp; '-----Special-----'];
for i=1:numel(optfields)
    temp=[temp; optfields{i}];
end
set(handles.plottype,'String',temp);
set(handles.plottype,'Value',1.0);
guiupdate(handles);

% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plottype.
function plottype_Callback(hObject, eventdata, handles)
% hObject    handle to plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plottype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plottype
contents = cellstr(get(hObject,'String'));
stemp=contents{get(hObject,'Value')};
set(handles.rtext,'String','Flux Surf');
set(handles.thetext,'String','Theta');
set(handles.tortext,'String','Zeta');
legend off;
if ~strcmp(stemp(1:1),'-')
    switch stemp
        case {'Summary'}
            handles.cuttype='text';
            set(handles.rcut,'Enable','off');
            set(handles.thetacut,'Enable','off');
            set(handles.zetacut,'Enable','off');
            set(handles.fluxcut,'Enable','off');
            set(handles.polcut,'Enable','off');
            set(handles.rzcut,'Enable','off');
            set(handles.threedcut,'Enable','off');
        case {'iotaf','mass','presf','beta_vol','phip','buco','bvco',...
                'phi','vp','overr','jcuru','jcurv',...
                'Jcuru','Jcurv','specw','Dmerc','Dshear','Dwell',...
                'Dcurr','Dgeod','jdotb','bdotgradv',...
                'dpdr','resid','nisla','correc','islw',...
                'correc2','edgeo','jdev','|B| Modes (axi)',...
                '|B| Modes (nonaxi)','pres_step','iota_step' 'fits',...
                'errp','pmap','omega','tpotb'}
            % Set Cuttype
            handles.cuttype='1d';
            % Turn off Cutplane options
            set(handles.cutplane,'SelectedObject',handles.rcut);
            set(handles.rcut,'Enable','on');
            set(handles.thetacut,'Enable','off');
            set(handles.zetacut,'Enable','off');
            set(handles.fluxcut,'Enable','off');
            set(handles.polcut,'Enable','off');
            set(handles.rzcut,'Enable','off');
            set(handles.threedcut,'Enable','off');
        case {'Flux Surface'}
            handles.cuttype='other';
            set(handles.cutplane,'SelectedObject',handles.rzcut);
            set(handles.rcut,'Enable','off');
            set(handles.thetacut,'Enable','off');
            set(handles.zetacut,'Enable','off');
            set(handles.fluxcut,'Enable','off');
            set(handles.polcut,'Enable','off');
            set(handles.rzcut,'Enable','on');
            set(handles.threedcut,'Enable','on');
        case 'Field Lines'
            handles.cuttype='other';
            set(handles.cutplane,'SelectedObject',handles.threedcut);
            set(handles.rcut,'Enable','off');
            set(handles.thetacut,'Enable','off');
            set(handles.zetacut,'Enable','off');
            set(handles.fluxcut,'Enable','off');
            set(handles.polcut,'Enable','off');
            set(handles.rzcut,'Enable','off');
            set(handles.threedcut,'Enable','on');
            set(handles.thetext,'String','nlines');
        case {'B-Field','J-Field','Flux Web','B-Field(pol)',...
                'Poincare (real)','Poincare (polar)'}
            handles.cuttype='zeta2';
            set(handles.rcut,'Enable','off');
            set(handles.thetacut,'Enable','off');
            set(handles.zetacut,'Enable','off');
            set(handles.fluxcut,'Enable','off');
            set(handles.polcut,'Enable','off');
            set(handles.rzcut,'Enable','on');
            set(handles.threedcut,'Enable','off');
            set(handles.rslide,'Enable','off');
            set(handles.theslide,'Enable','off');
            set(handles.torslide,'Enable','on');
            set(handles.cutplane,'SelectedObject',handles.rzcut);
        case 'LPK Plot'
            handles.cuttype='zeta2';
            set(handles.rcut,'Enable','off');
            set(handles.thetacut,'Enable','off');
            set(handles.zetacut,'Enable','off');
            set(handles.fluxcut,'Enable','off');
            set(handles.polcut,'Enable','off');
            set(handles.rzcut,'Enable','on');
            set(handles.threedcut,'Enable','off');
            set(handles.rslide,'Enable','off');
            set(handles.theslide,'Enable','off');
            set(handles.torslide,'Enable','off');
            set(handles.cutplane,'SelectedObject',handles.rzcut);
        case {'rbc' 'rbs' 'zbc' 'zbs'}
            handles.cuttype='r2';
            set(handles.rcut,'Enable','off');
            set(handles.thetacut,'Enable','off');
            set(handles.zetacut,'Enable','off');
            set(handles.fluxcut,'Enable','on');
            set(handles.polcut,'Enable','off');
            set(handles.rzcut,'Enable','off');
            set(handles.threedcut,'Enable','off');
            set(handles.rslide,'Enable','on');
            set(handles.theslide,'Enable','off');
            set(handles.torslide,'Enable','off');
            set(handles.cutplane,'SelectedObject',handles.fluxcut);
        otherwise
            switch handles.cuttype
                case {'1d','text'}
                    handles.cuttype='r1d';
                    set(handles.rcut,'Value',1.0);
                set(handles.cutplane,'SelectedObject',handles.rcut);
                case 'other'
                    handles.cuttype='zeta2';
                    set(handles.rzcut,'Value',1.0);
                    set(handles.cutplane,'SelectedObject',handles.rzcut);
            end
            set(handles.rcut,'Enable','on');
            set(handles.thetacut,'Enable','on');
            set(handles.zetacut,'Enable','on');
            set(handles.fluxcut,'Enable','on');
            set(handles.polcut,'Enable','on');
            set(handles.rzcut,'Enable','on');
            set(handles.threedcut,'Enable','on');
    end
    % Now we handle the zeta slider
    if handles.data.ntor == 0
        set(handles.torslide,'Enable','off');
    end
    guidata(hObject,handles);
    guiupdate(handles);
end

% --- Executes during object creation, after setting all properties.
function plottype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function torslide_Callback(hObject, eventdata, handles)
% hObject    handle to torslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.zetaval=int32(get(hObject,'Value'));
guidata(hObject,handles);
update_plots(handles);

% --- Executes during object creation, after setting all properties.
function torslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to torslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function update_plots(handles)
% UPDATE_PLOTS(handles) Updates the axes object
% This function updates the axes object
cla;
% Handle theta now being 3D
temp_theta=handles.theta;
if (length(size(handles.theta))==3)
    handles.theta=squeeze(handles.theta(handles.rval,:,handles.zetaval));
end
contents = cellstr(get(handles.plottype,'String'));
switch contents{get(handles.plottype,'Value')}
% Test Plot
    case 'Summary'
        c1=5; c2=25; c3=55; c4=75;
        axis([0 100 0 100]);
        % Basic Run info
        text(c1,95,'Run: ');
        if isfield(handles.data,'ns')
            text(c1,90,'NS:');
            text(c2,90,num2str(handles.data.ns));
        end
        if isfield(handles.data,'mpol')
            text(c1,85,'MPOL:');
            text(c2,85,num2str(handles.data.mpol));
        end
        if isfield(handles.data,'ntor')
            text(c1,80,'NTOR:');
            text(c2,80,num2str(handles.data.ntor));
        end
        if isfield(handles.data,'nfp')
            text(c1,75,'NFP:');
            text(c2,75,num2str(handles.data.nfp));
        end
        if isfield(handles.data,'input_extension')
            text(c2,95,handles.data.input_extension,...
                'Interpreter','none');
        end
        if strcmp(handles.data.datatype,'boozer')
            text(c1,70,'R:');
            text(c2,70,['[ ' num2str(handles.data.rmin) ...
                ' , ' num2str(handles.data.rmax) ' ]']);
        end
        if isfield(handles.data,'aspect')
            text(c1,60,'Aspect:');
            text(c2,60,num2str(handles.data.aspect));
        end
        if isfield(handles.data,'Rmajor')
            text(c1,55,'Rmajor:');
            text(c2,55,num2str(handles.data.Rmajor));
        end
        if isfield(handles.data,'Aminor')
            text(c1,50,'Aminor:');
            text(c2,50,num2str(handles.data.Aminor));
        end
        if isfield(handles.data,'Volume') || isfield(handles.data,'volume')
            text(c1,45,'Volume:');
            text(c2,45,num2str(handles.data.Volume));
        end
        if isfield(handles.data,'mgrid_file')
            text(c1,70,'MGRID\_FILE:');
            text(c2,70,handles.data.mgrid_file,...
                'Interpreter','none');
        end
        if isfield(handles.data,'betatot') || isfield(handles.data,'beta')
            text(c1,35,'\beta_{total}:');
            text(c2,35,num2str(handles.data.betatot));
        end
            % Beta Info
        if isfield(handles.data,'betapol')
            text(c1,30,'\beta_{poloidal}:');
            text(c2,30,num2str(handles.data.betapol));
        end
        if isfield(handles.data,'betator')
            text(c1,25,'\beta_{toroidal}:');
            text(c2,25,num2str(handles.data.betator));
        end
        if isfield(handles.data,'betaxis')
            text(c1,20,'\beta_{axis}:');
            text(c2,20,num2str(handles.data.betaxis));
        end
            % Magnetic Info
        if isfield(handles.data,'VolAvgB')
            text(c3,90,'<B>:');
            text(c4,90,num2str(handles.data.VolAvgB));
        end
        if isfield(handles.data,'b0')
            text(c3,85,'B_0:');
            text(c4,85,num2str(handles.data.b0));
        end
        if isfield(handles.data,'Itor')
            text(c3,80,'I_{toroidal}:');
            text(c4,80,num2str(handles.data.Itor));
        end
        if isfield(handles.data,'IonLarmor')
            text(c3,75,'Ion Larmor:');
            text(c4,75,num2str(handles.data.IonLarmor));
        end
        xlabel('');
        ylabel('');
        title('Summary of Run');
% 1D Arrays
    case 'iotaf'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.iotaf,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.iotaf,'k')
        xlabel(handles.xlabel);
        ylabel('Iota');
        title('Rotational Transform');
        % Code to overplot resonances
        if isfield(handles.data,'s_res')
            y_lims=[min(handles.data.iotaf) max(handles.data.iotaf)];
            hold on
            for i=1:size(handles.data.s_res,1)
                for j=1:size(handles.data.s_res,2);
                    if (handles.data.s_res(i,j) > 0)
                        plot(handles.data.s_res(i,j).*[1 1]./handles.data.ns,y_lims,'r');
                        text(handles.data.s_res(i,j)./handles.data.ns,y_lims(2),...
                            [num2str(handles.data.N(i),'%2d') '/' ...
                            num2str(handles.data.M(j),'%2d')]);
                    end
                end
            end
            hold off
        end
    case 'presf'
        if isfield(handles.data,'pres_plot')
            plot(handles.data.phi_plot,handles.data.pres_plot,'k');
            ylim([0 max(handles.data.pres_plot)]);
        else
            plot(handles.data.phi./handles.data.phi(handles.data.ns),...
                handles.data.presf,'k');
        end
        %plot(0:1./(handles.data.ns-1):1,handles.data.presf,'k')
        xlabel(handles.xlabel);
        ylabel('Pressure [Pa]');
        title('Pressure');
    case 'iota'
        plot(0:1./(handles.data.ns-1):1,handles.data.iota,'k')
        xlabel(handles.xlabel);
        ylabel('Iota');
        title('Rotational Transform');
        % Code to overplot resonances
        if isfield(handles.data,'s_res')
            y_lims=[min(handles.data.iota) max(handles.data.iota)];
            hold on
            for i=1:size(handles.data.s_res,1)
                for j=1:size(handles.data.s_res,2);
                    if (handles.data.s_res(i,j) > 0)
                        plot(handles.data.s_res(i,j).*[1 1]./handles.data.ns,y_lims,'r');
                        text(handles.data.s_res(i,j)./handles.data.ns,y_lims(2),...
                            [num2str(handles.data.N(i),'%2d') '/' ...
                            num2str(handles.data.M(j),'%2d')]);
                    end
                end
            end
            hold off
        end
    case 'pres'
        plot(0:1./(handles.data.ns-1):1,handles.data.pres,'k')
        xlabel(handles.xlabel);
        ylabel('Pressure [Pa]');
        title('Pressure');
    case 'beta_vol'
        plot(0:1./(handles.data.ns-1):1,handles.data.beta_vol,'k')
        xlabel(handles.xlabel);
        ylabel('Beta');
        title('Plasma Beta');
    case 'phip'
        if isfield(handles.data,'phipf')
            plot(0:1./(handles.data.ns-1):1,handles.data.phipf,'k')
        else
            plot(0:1./(handles.data.ns-1):1,handles.data.phip,'k')
        end
        xlabel(handles.xlabel);
        ylabel('Ploidal Flux [Wb]');
        title('Phip');
    case 'buco'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.buco,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.buco,'k')
        xlabel(handles.xlabel);
        ylabel('Buco');
        title('Buco');
    case 'bvco'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.bvco,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.bvco,'k')
        xlabel(handles.xlabel);
        ylabel('Bvco');
        title('Bvco');
    case 'phi'
        plot(0:1./(handles.data.ns-1):1,handles.data.phi,'k')
        xlabel(handles.xlabel);
        ylabel('Phi');
        title('Phi');
    case 'vp'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.vp,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.vp,'k')
        xlabel(handles.xlabel);
        ylabel('Vp');
        title('Vp');
    case 'overr'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.overr,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.overr,'k')
        xlabel(handles.xlabel);
        ylabel('Overr');
        title('Overr');
    case 'jcuru'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.jcuru,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.jcuru,'k')
        xlabel(handles.xlabel);
        ylabel('Jcuru');
        title('Jcuru');
    case 'jcurv'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.jcurv,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.jcurv,'k')
        xlabel(handles.xlabel);
        ylabel('Jcurv');
        title('Jcurv');
    case 'specw'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.specw,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.specw,'k')
        xlabel(handles.xlabel);
        ylabel('Specw');
        title('Specw');
    case 'Dmerc'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.Dmerc,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.Dmerc,'k')
        xlabel(handles.xlabel);
        ylabel('DMerc');
        title('DMerc');
    case 'Dshear'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.Dshear,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.Dshear,'k')
        xlabel(handles.xlabel);
        ylabel('D \iota / D \phi');
        title('Dshear');
    case 'Dwell'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.Dwell,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.Dwell,'k')
        xlabel(handles.xlabel);
        ylabel('Dwell');
        title('Dwell');
    case 'Dcurr'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.Dcurr,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.Dcurr,'k')
        xlabel(handles.xlabel);
        ylabel('Dcurr');
        title('Dcurr');
    case 'Dgeod'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.Dgeod,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.Dgeod,'k')
        xlabel(handles.xlabel);
        ylabel('Dgeod');
        title('Dgeod');
    case 'jdotb'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.jdotb,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.jdotb,'k')
        xlabel(handles.xlabel);
        ylabel('J \cdot B');
        title('J \cdot B');
    case 'bdotgradv'
        plot(handles.data.phi./handles.data.phi(handles.data.ns),...
            handles.data.bdotgradv,'k');
        %plot(0:1./(handles.data.ns-1):1,handles.data.bdotgradv,'k')
        xlabel(handles.xlabel);
        ylabel('B \cdot \nabla v');
        title('B \cdot \nabla v');
    case 'dpdr'
        plot(0:1./(handles.data.ns-1):1,handles.data.dpdr,'k')
        xlabel(handles.xlabel);
        ylabel('Dp/Dr');
        title('Radial Pressure Derivative');
    case 'beta'
        plot(0:1./(handles.data.ns-1):1,handles.data.beta,'k')
        xlabel(handles.xlabel);
        ylabel('<\beta>');
        title('Plasma Beta');
    case 'resid'
        plot(0:handles.data.niter,handles.data.resid,'k')
        xlabel('Iteration');
        ylabel('Force Residual');
        title('Normalized Force Residuals');
    case 'correc'
        plot(0:handles.data.niter,handles.data.correc,'k')
        set(gca,'YScale','log')
        xlabel('Iteration');
        ylabel('Correction');
        title('Correction (correc)');
    case 'correc2'
        plot(0:handles.data.niter,handles.data.correc2,'k')
        set(gca,'YScale','log')
        xlabel('Iteration');
        ylabel('Correction^2');
        title('Correction^2 (correc2)');
    case 'errp'
        nrho=size(handles.data.errp,1);
        plot(handles.data.itermax.*(0:1/(nrho-1):1),handles.data.errp(:,1:3))
        set(gca,'YScale','log')
        xlabel('Iteration');
        ylabel('Error');
        title('Error vs. Iteration');
        legend('R','PSI','AX');
    case 'islw'
        plot(0:handles.data.niter,handles.data.islw,'k')
        xlabel('Iteration');
        ylabel('Width');
        title('Island Width');
    case 'jdev'
        bar(1:handles.data.k,handles.data.jdev,'k')
        rhodex=0:handles.data.k;
        xlabel({'Radial Surface '; num2str(rhodex(handles.data.jdev==1))});
        ylabel('Bad Surfaces');
        title('Bad Surfaces');
    case 'nisla'
        bar(0:handles.data.niter,handles.data.nisla,'k')
        hold on
        bar(0:handles.data.niter,handles.data.nislr,'r','BarWidth',0.6);
        bar(0:handles.data.niter,handles.data.nislha,'b','BarWidth',0.4);
        xlabel('Iteration');
        ylabel('Number of Islands');
        title('Island Inventory');
        legend('Total Islands','Resolved Islands','Hudson Islands');
        hold off
    case 'edgeo'
        inner=handles.data.edgei(handles.data.edgei > 0);
        outer=handles.data.edgeo(handles.data.edgeo > 0);
        center=(inner+outer)./2;
        errorbar(1:length(center),center,...
            abs(inner-center),abs(outer-center),'.k');
        xlabel('Island Number');
        ylabel('Radial extent');
        ylim([0 1]);
        title('Island Edges');
    case 'mu'
        plot(handles.data.phi_plot,handles.data.mu_plot,'k');
        xlabel('Normalized Toroidal Flux');
        ylabel('Helicity Multiplier');
        title('Stepped Helicity Multiplier Profile');
    case 'iota_step'
        phi = [0; handles.data.phi_step];
        bar(phi,[0; handles.data.iota_step],0.1);
        xlabel('Normalized Toroidal Flux');
        ylabel('Rotational Transform');
        title('Interface Rotational Transform');
    case 'fits'
        plot(handles.data.fits(1,:),log(handles.data.fits(3,:)));
        xlabel('Iteration');
        ylabel('log(Fbal)');
        title('Force Balance Evolution');
    case 'pmap'
        plot(0:1./(handles.data.ns-1):1,handles.data.pmap,'k')
        xlabel(handles.xlabel);
        ylabel('<p>');
        title('Flux Surface Averaged Pressure');
    case 'omega'
        plot(0:1./(handles.data.ns-1):1,handles.data.omega,'k')
        xlabel(handles.xlabel);
        ylabel('Freq.');
        title('Toroidal Angular Frequency');
    case 'tpotb'
        plot(0:1./(handles.data.ns-1):1,handles.data.tpotb,'k')
        xlabel(handles.xlabel);
        ylabel('Temperature');
        title('Normalized Temperature');
        
% 3D Arrays
    case {'b','lam','p','g','bs','bu','bv','b_s','b_u','b_v','br','bphi','bz',...
            'currs' 'curru','currv','jr','jphi','jz','brho','btheta','bnorm',...
            'jpara','press','fnx','fny','fnz','fn','jxbx','jxby','jxbz',...
            'dpds','ptran','prot','prpr','U_rot','vphi'}
        % Get Title Strings
        string=contents{get(handles.plottype,'Value')};
        if strcmp(string,'b'),name='|B|';
        elseif strcmp(string,'g'),name='g';
        elseif strcmp(string,'bs'),name='B^s';
        elseif strcmp(string,'bu'),name='B^u';
        elseif strcmp(string,'bv'),name='B^v';
        elseif strcmp(string,'b_s'),name='B_s';
        elseif strcmp(string,'b_u'),name='B_u';
        elseif strcmp(string,'b_v'),name='B_v';
        elseif strcmp(string,'br'),name='B_R';
        elseif strcmp(string,'bphi'),name='B_\phi';
        elseif strcmp(string,'bz'),name='B_Z';
        elseif strcmp(string,'curru'),name='j^u';
        elseif strcmp(string,'currv'),name='j^v';
        elseif strcmp(string,'jr'),name='j_r';
        elseif strcmp(string,'jphi'),name='j_\phi';
        elseif strcmp(string,'jz'),name='j_z';
        elseif strcmp(string,'brho'),name='B_\rho';
        elseif strcmp(string,'btheta'),name='B_\theta';
        elseif strcmp(string,'bnorm'),name='Normal B-Field';
        elseif strcmp(string,'ptran'),name='Boozer p transform';
        elseif strcmp(string,'prot'),name='Pressure';
        elseif strcmp(string,'U_rot'),name='Rotational Energy';
        elseif strcmp(string,'vphi'),name='Toroidal Velocity';
        elseif strcmp(string,'jpara')
            if isfield(handles.data,'jlc_title')
                name=handles.data.jlc_title;
            else
                name='J Parallel';
            end
        elseif strcmp(string,'press')
            if isfield(handles.data,'press_title')
                name=handles.data.press_title;
            else
                name='Plasma Pressure';
            end
        else name=string;
        end
        % Now Plot
        f=handles.(string);
        if strcmp(handles.cuttype,'r1d')
            plot(f(:,handles.thetaval,handles.zetaval))
            xlabel(handles.xlabel);
            ylabel(name);
            title([name ' at theta=',...
                num2str(handles.theta(handles.thetaval)),...
                ' \zeta=',...
                num2str(handles.zeta(handles.zetaval))]);
            % Code to overplot resonances
            y_lims=[min(f(:,handles.thetaval,handles.zetaval)) max(f(:,handles.thetaval,handles.zetaval))];
            if isfield(handles.data,'s_res')
                hold on
                for i=1:size(handles.data.s_res,1)
                    for j=1:size(handles.data.s_res,2);
                        if (handles.data.s_res(i,j) > 0)
                            plot(handles.data.s_res(i,j).*[1 1],y_lims,'r');
                            text(handles.data.s_res(i,j),y_lims(2),...
                                [num2str(handles.data.N(i),'%2d') '/' ...
                                num2str(handles.data.M(j),'%2d')]);
                        end
                    end
                end
                hold off
            end
        elseif strcmp(handles.cuttype,'theta1d')
            plot(handles.theta,...
                f(handles.rval,:,handles.zetaval))
            xlabel('Theta (\theta) [rad]');
            ylabel(name);
            title([name 'on flux surface r=',num2str(handles.rval),...
                ' at \zeta=',num2str(handles.zeta(handles.zetaval))]);
        elseif strcmp(handles.cuttype,'zeta1d')
            plot(squeeze(handles.zeta),...
                squeeze(f(handles.rval,handles.thetaval,:)))
            xlabel('Zeta (\zeta) [rad]');
            ylabel(name);
            title([name ' on flux surface r=',num2str(handles.rval),...
                ' at theta=',num2str(handles.theta(handles.thetaval))]);
        elseif strcmp(handles.cuttype,'r2')
            if handles.data.nfp==1
                h=pcolor(handles.zeta,handles.theta,squeeze(f(handles.rval,:,:)));
            else
                zetafp=handles.zeta;
                h=pcolor(zetafp,handles.theta,...
                    squeeze(f(handles.rval,:,:)));
            end
            set(h,'EdgeColor','none');
            xlabel('Zeta (\zeta) [rad]');
            ylabel('Theta (\theta) [rad]');
            title([name ' on flux surface=',...
                num2str(handles.rval)]);
            % Overplot fieldlines if Boozer Coordinates
            if strcmp(handles.data.datatype,'boozer') && strcmp(name,'|B|') && (handles.data.iota(handles.rval) > 0)
                hold on
                    if handles.data.nfp==1
                        plot(handles.zeta,handles.zeta.*handles.data.iota(handles.rval),'k');
                    else
                        dtheta=max(zetafp.*handles.data.iota(handles.rval));
                        nlines=round(2*pi/dtheta);
                        for i=1:nlines
                            fieldline=dtheta*(i-1)+zetafp.*handles.data.iota(handles.rval);
                            dex=(fieldline<2*pi);
                            plot(zetafp(dex),fieldline(dex),'k');
                        end
                    end
                hold off
            end
        elseif strcmp(handles.cuttype,'theta2')
            if handles.data.nfp==1
                h=pcolor(repmat(handles.zeta,[handles.data.ns 1]),...
                    1:handles.data.ns,...
                    squeeze(f(:,handles.thetaval,:)));
            else
                zetafp=handles.zeta;
                h=pcolor(repmat(zetafp,[handles.data.ns 1]),...
                    1:handles.data.ns,...
                    squeeze(f(:,handles.thetaval,:)));
            end
            set(h,'EdgeColor','none');
            xlabel('Zeta (\zeta)');
            ylabel('Radial Grid');
            title([name ' for \theta=',...
                num2str(handles.theta(handles.thetaval))]);
            axis tight
        elseif strcmp(handles.cuttype,'zeta2')
            torocont(handles.r,handles.z,f,handles.zetaval)
            xlabel('Radius (R) [m]');
            ylabel('Elevation (Z) [m]');
            title([name ' at \phi=',num2str(handles.zeta(handles.zetaval))]);
        elseif strcmp(handles.cuttype,'3D')
            set(handles.rtext,'String','Flux');
            r = repmat(handles.r(:,:,1:end-1),[1 1 handles.data.nfp]);
            z = repmat(handles.z(:,:,1:end-1),[1 1 handles.data.nfp]);
            f = repmat(f(:,:,1:end-1),[1 1 handles.data.nfp]);
            r(:,:,end+1) = r(:,:,1);
            z(:,:,end+1) = z(:,:,1);
            f(:,:,end+1) = f(:,:,1);
            phi = 0:2*pi/(size(r,3)-1):2*pi;
            if (handles.rval > 1)
                isotoro(r,z,phi,handles.rval,f);
                title([name ' on Flux Surface (ns=',num2str(handles.rval),')']);
            else
                plot3(squeeze(r(1,1,:)).*cos(phi'),...
                    squeeze(r(1,1,:)).*sin(phi'),...
                    squeeze(z(1,1,:)),'k')
                title('Magnetic Axis');
            end
            xlabel('X [m]');
            ylabel('Y [m]');
            zlabel('Z [m]');
        end
    % Fourier Arrays
    case {'rbc' 'rbs' 'zbc' 'zbs'}
        % Get Title Strings
        string=contents{get(handles.plottype,'Value')};
        if strcmp(string,'rbc'),name='EVEN R Coefficients (cos)';
        elseif strcmp(string,'rbs'),name='ODD R Coefficients (sin)';
        elseif strcmp(string,'zbc'),name='EVEN Z Coefficients (cos)';
        elseif strcmp(string,'zbs'),name='ODD Z Coefficients (sin)';
        else name=string;
        end
        % Now Plot
        f=handles.(string);
        if ~(handles.data.nfp==1)
            surf(0:handles.data.mpol-1,-handles.data.ntor:handles.data.ntor,f(handles.rval,:,:).^2)
            xlabel('Poloidal Mode Number (m)');
            ylabel('Toroidal Mode Number (n)');
            zlabel('Fourier Amplitude');
            set(gca,'ZScale','log');
        else
            plot(0:handles.data.mpol-1,f(handles.rval,:).^2);
            set(gca,'YScale','log');
            ylabel('Fourier Amplitude');
            xlabel('Poloidal Mode Number (m)');
        end
        title([name ' (ns=',num2str(handles.rval),')']);
% Special Plots
    case 'Flux Web'
        colorbar('off');
        if get(handles.rzcut,'Value')==1
            set(handles.rslide,'Enable','off');
            set(handles.theslide,'Enable','off');
            set(handles.torslide,'Enable','on');
            handles.cuttype='zeta2';
            surfs=2:handles.data.ns;
        else
            handles.cuttype='3D';
            set(handles.rslide,'Enable','on');
            set(handles.theslide,'Enable','off');
            set(handles.torslide,'Enable','off');
        end
        if strcmp(handles.cuttype,'zeta2')
            hp=pcolor(handles.r(:,:,handles.zetaval),...
                handles.z(:,:,handles.zetaval),...
                0.0*handles.r(:,:,handles.zetaval));
            set(hp,'FaceColor','none','EdgeColor','black');
            hold on
            plot(handles.r(1,1,handles.zetaval),handles.z(1,1,handles.zetaval),'o');
            hold off
            xlabel('Radius (R) [m]');
            ylabel('Elevation (Z) [m]');
            title(strcat('Flux Surfaces at zeta=',num2str(handles.zeta(handles.zetaval))));
            axis equal
        elseif strcmp(handles.cuttype,'3D')
            set(handles.rtext,'String','Flux');
            isotoro(handles.r,handles.z,handles.zeta,handles.rval);
            title(strcat('Flux Surface (ns=',num2str(handles.rval),')'));
            xlabel('X [m]');
            ylabel('Y [m]');
            zlabel('Z [m]');
        end
    case 'B_R'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','on');
        pixplot(handles.data.raxis,handles.data.zaxis,squeeze(handles.data.B_R(:,handles.zetaval,:)))
        xlabel('R [m]');
        ylabel('Z [m]');
        handles.cuttype='zeta2';
    case 'B_PHI'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','on');
        pixplot(handles.data.raxis,handles.data.zaxis,squeeze(handles.data.B_PHI(:,handles.zetaval,:)))
        xlabel('R [m]');
        ylabel('Z [m]');
        handles.cuttype='zeta2';
    case 'B_Z'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','on');
        pixplot(handles.data.raxis,handles.data.zaxis,squeeze(handles.data.B_Z(:,handles.zetaval,:)))
        xlabel('R [m]');
        ylabel('Z [m]');
        handles.cuttype='zeta2';
    case 'S_ARR'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','on');
        pixplot(handles.data.raxis,handles.data.zaxis,squeeze(handles.data.S_ARR(:,handles.zetaval,:)))
        xlabel('R [m]');
        ylabel('Z [m]');
        handles.cuttype='zeta2';
    case 'Flux Surface'
        colorbar('off');
        if get(handles.rzcut,'Value')==1
            set(handles.rslide,'Enable','on');
            set(handles.theslide,'Enable','off');
            set(handles.torslide,'Enable','on');
            handles.cuttype='zeta2';
            if handles.rval==1
                surfs=[2 handles.data.ns];
            else
                surfs=2:(handles.data.ns-2)/double(handles.rval-1):handles.data.ns;
                surfs=round(surfs);
            end
        else
            handles.cuttype='3D';
            set(handles.rslide,'Enable','on');
            set(handles.theslide,'Enable','off');
            set(handles.torslide,'Enable','off');
        end
        if strcmp(handles.cuttype,'zeta2')
            set(handles.rtext,'String','ns');
            set(handles.rslide,'Enable','on');
            toroslice(handles.r,handles.zetaval,handles.z,surfs);
            hold on
            plot(handles.r(1,1,handles.zetaval),handles.z(1,1,handles.zetaval),'+');
            if (isfield(handles.data,'ideal_arr'))
                for ik=1:length(handles.data.ideal_arr)
                    s_temp=handles.data.ideal_arr(ik);
                    plot(handles.r(s_temp,:,handles.zetaval),handles.z(s_temp,:,handles.zetaval),'LineWidth',2.0,'Color','k');
                end
            end
            hold off
            xlabel('Radius (R) [m]');
            ylabel('Elevation (Z) [m]');
            title(strcat('Flux Surfaces at zeta=',num2str(handles.zeta(handles.zetaval))));
            axis equal
        elseif strcmp(handles.cuttype,'3D')
            
            set(handles.rtext,'String','Flux');
            r = repmat(handles.r(:,:,1:end-1),[1 1 handles.data.nfp]);
            z = repmat(handles.z(:,:,1:end-1),[1 1 handles.data.nfp]);
            r(:,:,end+1) = r(:,:,1);
            z(:,:,end+1) = z(:,:,1);
            phi = 0:2*pi/(size(r,3)-1):2*pi;
            if (handles.rval > 1)
                isotoro(r,z,phi,handles.rval);
                title(strcat('Flux Surface (ns=',num2str(handles.rval),')'));
            else
                plot3(squeeze(r(1,1,:)).*cos(phi'),...
                    squeeze(r(1,1,:)).*sin(phi'),...
                    squeeze(z(1,1,:)),'k')
                title('Magnetic Axis');
            end
            xlabel('X [m]');
            ylabel('Y [m]');
            zlabel('Z [m]');
            camlight left;
        end
    case 'Field Lines' % doesn't work right in straight field line coords
        set(handles.rslide,'Enable','on');
        set(handles.theslide,'Enable','on');
        set(handles.torslide,'Enable','off');
        set(handles.thetext,'String','nlines');
        if handles.rval < 2
            handles.rval=2;
        end
        nt=1:round(size(handles.theta,2)/handles.thetaval):size(handles.theta,2);
        torolines_modb(handles.r(:,nt,:),handles.z(:,nt,:),-handles.zeta,[handles.rval],handles.b);
        title(strcat('Field Lines (on surface ns=',num2str(handles.rval),')'));
        xlabel('X [m]');
        ylabel('Y [m]');
        zlabel('Z [m]');
    case 'B-Field'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','on');
        torocont(handles.r,handles.z,handles.bphi,handles.zetaval);
        % Note we normalize the quiver plot to bphi
        hold on
        quiver(handles.r(:,:,handles.zetaval),...
            handles.z(:,:,handles.zetaval),...
            handles.br(:,:,handles.zetaval),...
            handles.bz(:,:,handles.zetaval),...
            'Color','black');
        hold off
        xlabel('R [m]');
        ylabel('Z [m]');
        title(strcat('Magnetic Field at zeta=',num2str(handles.zeta(handles.zetaval))));
        handles.cuttype='other';
    case 'B-Field(pol)'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','on');
        torocont(handles.r,handles.z,handles.bphi,handles.zetaval);
        % Note we normalize the quiver plot to bphi
        hold on
        quiver(handles.r(:,:,handles.zetaval),...
            handles.z(:,:,handles.zetaval),...
            handles.br_pol(:,:,handles.zetaval),...
            handles.bz_pol(:,:,handles.zetaval),...
            'Color','black');
        hold off
        xlabel('R [m]');
        ylabel('Z [m]');
        title(strcat('Poloidal Magnetic Field at zeta=',num2str(handles.zeta(handles.zetaval))));
        handles.cuttype='other';
    case 'J-Field'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','on');
        torocont(handles.r,handles.z,handles.jphi,handles.zetaval);
        % Note we normalize the quiver plot to bphi
        hold on
        quiver(handles.r(:,:,handles.zetaval),...
            handles.z(:,:,handles.zetaval),...
            handles.jr(:,:,handles.zetaval),...
            handles.jz(:,:,handles.zetaval),...
            'Color','black');
        hold off
        xlabel('R [m]');
        ylabel('Z [m]');
        title(strcat('Current Density at zeta=',num2str(handles.zeta(handles.zetaval))));
        handles.cuttype='other';
    case 'SPEC Grid'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','on');
        plot_spec_grid(handles.data,handles.zetaval);
        xlabel('R [m]');
        ylabel('Z [m]');
        title('SPEC Finite Element Grid');
        set(handles.torslide,'Max',handles.data.nzeta_grid-1);
        set(handles.torslide,'Min',1.0);
        set(handles.torslide,'SliderStep',[1.0 1.0]./(handles.data.nzeta_grid-1));
        handles.cuttype='zeta2';
    case 'Poincare (real)'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        if (size(handles.data.R_lines,2) == 1)
            set(handles.torslide,'Enable','off');
            handles.zetaval = 1;
        else
            set(handles.torslide,'Enable','on');
        end
        R=squeeze(handles.data.R_lines(:,handles.zetaval,:));
        Z=squeeze(handles.data.Z_lines(:,handles.zetaval,:));
        plot(R',Z','.');
        hold on;
        rho_dex=[1 handles.data.ideal_arr];
        plot(squeeze(handles.r(rho_dex,:,handles.zetaval))',squeeze(handles.z(rho_dex,:,handles.zetaval))','k');
        hold off;
        xlabel('R [m]');
        ylabel('Z [m]');
        title('Poincare Plot');
        handles.cuttype='zeta2';
        axis equal
        handles.ntor = handles.data.npoinc;
        %set(handles.torslide,'Max',handles.data.nzeta_grid);
        %set(handles.torslide,'Min',1.0);
        %step = 1.0/handles.data.nzeta_grid;
        %set(handles.torslide,'SliderStep',[step step]);
    case 'Poincare (polar)'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        if (length(size(handles.data.R_lines))==2)
            set(handles.torslide,'Enable','off');
            handles.zetaval = 1;
            rho=handles.data.rho_lines(:,:)';
            theta=handles.data.th_lines(:,:)';
        else
            set(handles.torslide,'Enable','on');
            rho=squeeze(handles.data.rho_lines(:,handles.zetaval,:))';
            theta=squeeze(handles.data.th_lines(:,handles.zetaval,:))';
        end
        plot(rho.*cos(theta),rho.*sin(theta),'.','MarkerSize',0.1)
        xlabel('R [m]');
        ylabel('Z [m]');
        title('Poincare Plot (polar)');
        handles.cuttype='poincare_polar';
        axis equal
    case 'Poincare (rho)'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        if (length(size(handles.data.R_lines))==2)
            set(handles.torslide,'Enable','off');
            handles.zetaval = 1;
            rho=handles.data.rho_lines(:,:)';
            theta=handles.data.th_lines(:,:)';
        else
            set(handles.torslide,'Enable','on');
            rho=squeeze(handles.data.rho_lines(:,handles.zetaval,:))';
            theta=squeeze(handles.data.th_lines(:,handles.zetaval,:))';
        end
        plot(theta,rho,'.','MarkerSize',0.5);
        %plot(rho.*cos(theta),rho.*sin(theta),'.')
        ylabel('rho [norm]');
        xlabel('theta [rad]');
        title('Poincare Plot (rho-theta)');
        handles.cuttype='poincare_polar';
        axis tight
    case 'LPK Plot'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','off');
        surfs=handles.data.ns;
        nzeta=size(handles.zeta,2);
        cuts=[1 round(nzeta/4)+1 round(nzeta/2)+1];
        toroslice(handles.r,cuts(1),handles.z,surfs,'b');
        hold on
        toroslice(handles.r,cuts(2),handles.z,surfs,'g');
        toroslice(handles.r,cuts(3),handles.z,surfs,'r');
        plot(squeeze(handles.r(1,1,cuts(1))),squeeze(handles.z(1,1,cuts(1))),'+b');
        plot(squeeze(handles.r(1,1,cuts(2))),squeeze(handles.z(1,1,cuts(2))),'+g');
        plot(squeeze(handles.r(1,1,cuts(3))),squeeze(handles.z(1,1,cuts(3))),'+r');
        hold off
        title('Flux Surface Evolution');
        xlabel('R [m]');
        ylabel('Z [m]');
        colorbar off
        handles.cuttype='other';
        axis equal
    case '|B| Modes (axi)'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','off');
        dex_n=handles.data.xn==0;
        dex_m=handles.data.xm==0;
        dex=logical(dex_n.*dex_m);
        plot(1:handles.data.ns,squeeze(handles.data.bmnc(dex,:)),'k');
        leg_text={'m = 0, n = 0'};
        hold on
        amplitudes=max(handles.data.bmnc,[],2);
        amplitudes(dex_m)=0;
        amplitudes(dex_n==0)=0;
        line_color=lines(4);
        for i=1:4
            max_amp=max(amplitudes);
            dex=find(amplitudes==max_amp,1,'first');
            leg_text=[leg_text; ['m = ' num2str(handles.data.xm(dex),'%2d') ', n = 0']];
            plot(1:handles.data.ns,handles.data.bmnc(dex,:),'Color',line_color(i,:));
            amplitudes(amplitudes >= max_amp) = 0;
        end
        title('Axisymmetric |B| modes');
        xlabel(handles.xlabel);
        ylabel('|B|');
        legend(leg_text,'Location','NorthWest');
        handles.cuttype='r1d';
        axis tight;
    case '|B| Modes (nonaxi)'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','off');
        dex_n=handles.data.xn==0;
        amplitudes=max(handles.data.bmnc,[],2);
        amplitudes(dex_n)=0;
        max_amp=max(amplitudes);
        dex=find(amplitudes==max_amp,1,'first');
        plot(1:handles.data.ns,handles.data.bmnc(dex,:),'k');
        leg_text={['m = ' num2str(handles.data.xm(dex),'%2d') ', n = ' num2str(handles.data.xn(dex))]};
        hold on
        line_color=lines(7);
        for i=1:7
            amplitudes(amplitudes>=max_amp)=0.0;
            max_amp=max(amplitudes);
            dex=find(amplitudes==max_amp,1,'first');
            plot(1:handles.data.ns,squeeze(handles.data.bmnc(dex,:)),'Color',line_color(i,:));
            leg_text=[leg_text; ['m = ' num2str(handles.data.xm(dex),'%2d') ', n = ' num2str(handles.data.xn(dex))]];
        end
        title('Non-axisymmetric |B| modes');
        xlabel(handles.xlabel);
        ylabel('|B|');
        legend(leg_text,'Location','NorthWest');
        handles.cuttype='r1d';
        axis tight;
    case '|B| Field Line'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','off');
        theta=handles.zeta'*handles.data.iota; %[zeta,ns]
        hold on;
        leg_text={};
        surfaces=[2:handles.data.ns/3:handles.data.ns handles.data.ns];
        line_color=lines(length(surfaces));
        for i=1:length(surfaces)
            t=surfaces(i);
            bline=zeros(1,size(theta,1));
            leg_text=[leg_text; ['Surface: ' num2str(t,'%3d')]];
            for n=1:size(theta,1)
                bline(n)=interp2(handles.theta,handles.zeta,squeeze(handles.b(t,:,:))',squeeze(theta(n,t)),handles.zeta(n));
            end
            plot(squeeze(theta(:,t)),bline,'Color',line_color(i,:));
        end
        title('|B| Along Field Lines');
        xlabel('Theta (\theta) [radians]');
        ylabel('|B|');
        legend(leg_text,'Location','NorthWest');
        handles.cuttype='r1d';
        axis tight;
        
end
switch handles.cuttype
    case 'text'
        zoom off
        colorbar('off')
        rotate3d off
    case {'1d','r1d','theta1d','zeta1d'}
        axis tight
        zoom on
        colorbar('off')
        rotate3d off
    case 'r2'
        axis tight
        zoom on
        colorbar('FontSize',16)
        rotate3d off
        colormap jet
    case 'theta2'
        axis tight
        zoom on
        colorbar
        rotate3d off
        colormap jet
    case {'zeta2','other'}
        set(handles.graph,'XLim',[handles.xmin handles.xmax]);
        set(handles.graph,'YLim',[handles.zmin handles.zmax]);
        zoom on
        colorbar('FontSize',16)
        rotate3d off
        if handles.plotves
            hold on
            ph=plot_vessel(handles.ves_data,'phi',handles.zeta(handles.zetaval));
            hold off
            set(ph,'Color','black');
        end
        if isfield(handles.data,'hitsrf')
            hold on
            ph=plot(handles.r(handles.data.hitsrf,:,handles.zetaval),...
                handles.z(handles.data.hitsrf,:,handles.zetaval),'--k');
            hold off
            switch contents{get(handles.plottype,'Value')}
                case{'Flux Surface','Flux Web'}
                    set(ph,'Color','red','LineStyle','-');
            end
        end 
        colormap jet
    case '3D'
        axis equal
        rotate3d on
        zoom off
        if handles.plotves
            hold on
            ph=plot_vessel(handles.ves_data,'wire');
            hold off
        end
        colormap jet
    case {'poincare','poincare_polar'}
        %axis equal
        zoom on
        %colorbar('FontSize',16)
        rotate3d off
        colormap jet
        
end
% Now clean up 'Special' Plots
switch contents{get(handles.plottype,'Value')}
    case {'Flux Surface','Poincare (real)','Poincare (polar)'}
        colorbar off
    case {'B-Field','J-Field','B-Field(pol)'}
        zoom on
        rotate3d off
        colorbar('FontSize',16)
    case 'Field Lines'
        zoom off
        rotate3d on
    case {'LPK Plot','Flux Web','|B| Modes (axi)'}
        zoom on
        rotate3d off
        colorbar off
end
% Now Fix theta
handles.theta=temp_theta;
        
        
function guiupdate(handles)
% GUIUPDATE(handles) Updates the GUI
% This function updates the GUI
switch handles.cuttype
    case 'text'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','off');
    case '1d'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','off');
    case 'r1d'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','on');
        set(handles.torslide,'Enable','on');
    case 'theta1d'
        set(handles.rslide,'Enable','on');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','on');
    case 'zeta1d'
        set(handles.rslide,'Enable','on');
        set(handles.theslide,'Enable','on');
        set(handles.torslide,'Enable','off');
    case 'r2'
        set(handles.rslide,'Enable','on');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','off');
    case 'theta2'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','on');
        set(handles.torslide,'Enable','off');
    case 'zeta2'
        set(handles.rslide,'Enable','off');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','on');
    case '3D'
        set(handles.rslide,'Enable','on');
        set(handles.theslide,'Enable','off');
        set(handles.torslide,'Enable','off');
end
if handles.data.ntor == 0
    set(handles.torslide,'Enable','off');
    set(handles.thetacut,'Enable','on');
    set(handles.zetacut,'Enable','off');
    set(handles.fluxcut,'Enable','off');
    set(handles.polcut,'Enable','off');
end
set(handles.ntheta,'String',num2str(handles.mpol));
set(handles.nzeta,'String',num2str(handles.ntor));
update_plots(handles);



function ntheta_Callback(hObject, eventdata, handles)
% hObject    handle to ntheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ntheta as text
%        str2double(get(hObject,'String')) returns contents of ntheta as a double
handles.mpol=str2double(get(hObject,'String'));
handles.theta=0:2*pi/(handles.mpol-1):2*pi;
set(handles.theslide,'Max',handles.mpol);
set(handles.theslide,'Min',1.0);
set(handles.theslide,'SliderStep',[1.0 1.0]./(handles.mpol-1));
set(handles.theslide,'Value',1.0);
guidata(hObject, handles);
transf(hObject,handles);
guiupdate(handles);

% --- Executes during object creation, after setting all properties.
function ntheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ntheta (see GCBO)
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
%        str2double(get(hObject,'String')) returns contents of nzeta as
%        a double
handles.ntor=str2double(get(hObject,'String'));
handles.zeta=0:2*pi/double(handles.ntor-1):2*pi;
set(handles.torslide,'Max',handles.ntor);
set(handles.torslide,'Min',1.0);
set(handles.torslide,'SliderStep',[1.0 1.0]./double(int32(handles.ntor)-1));
set(handles.torslide,'Value',1.0);
guidata(hObject, handles);
transf(hObject,handles);
guiupdate(handles);

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

function handles=transf(hObject,handles)
% Transforms Fourier Spectrums
cla;
axis([0 1 0 1]);
text(0.35,0.5,'Computing Values');
text(0.35,0.4,'Please Wait');
switch handles.data.datatype
    case 'wout'
        % Label X-axis
        handles.xlabel='Normalize Toroidal Flux (s)';
        % Reset theta to 1D
        handles.theta=0:2*pi/(handles.mpol-1):2*pi;
        set(handles.statustext,'String','Computing lambda');
        pause(.01);
        % Switch to field period representation
        xn = handles.data.xn./handles.data.nfp;
        xn_nyq = handles.data.xn_nyq./handles.data.nfp;
        % Transform to straight field line coords
        handles.l=sfunct(handles.theta,handles.zeta,handles.data.lmns,handles.data.xm,xn);
        if handles.data.iasym==1
            set(handles.statustext,'String','Asymetric VMEC detected!');
            pause(.01);
            set(handles.statustext,'String','Adding lbc to lambda');
            pause(.01);
            handles.l=handles.l-cfunct(handles.theta,handles.zeta,handles.data.lmnc,handles.data.xm,xn);
        end
        % Transform quantities
        set(handles.statustext,'String','Computing r');
        pause(.01);
        handles.r=cfunct(handles.theta,handles.zeta,handles.data.rmnc,handles.data.xm,xn);
        set(handles.statustext,'String','Computing z');
        pause(.01);
        handles.z=sfunct(handles.theta,handles.zeta,handles.data.zmns,handles.data.xm,xn);
        set(handles.statustext,'String','Computing lam');
        pause(.01);
        handles.lam=sfunct(handles.theta,handles.zeta,handles.data.lmns,handles.data.xm,xn);
        set(handles.statustext,'String','Computing B');
        pause(.01);
        handles.b=cfunct(handles.theta,handles.zeta,handles.data.bmnc,handles.data.xm_nyq,xn_nyq);
        set(handles.statustext,'String','Computing g');
        pause(.01);
        handles.g=cfunct(handles.theta,handles.zeta,handles.data.gmnc,handles.data.xm_nyq,xn_nyq);
        set(handles.statustext,'String','Computing b_s');
        pause(.01);
        handles.b_s=sfunct(handles.theta,handles.zeta,handles.data.bsubsmns,handles.data.xm_nyq,xn_nyq);
        set(handles.statustext,'String','Computing b_u');
        pause(.01);
        handles.b_u=cfunct(handles.theta,handles.zeta,handles.data.bsubumnc,handles.data.xm_nyq,xn_nyq);
        set(handles.statustext,'String','Computing b_v');
        pause(.01);
        handles.b_v=cfunct(handles.theta,handles.zeta,handles.data.bsubvmnc,handles.data.xm_nyq,xn_nyq);
        set(handles.statustext,'String','Computing b^u');
        pause(.01);
        handles.bu=cfunct(handles.theta,handles.zeta,handles.data.bsupumnc,handles.data.xm_nyq,xn_nyq);
        set(handles.statustext,'String','Computing b^v');
        pause(.01);
        handles.bv=cfunct(handles.theta,handles.zeta,handles.data.bsupvmnc,handles.data.xm_nyq,xn_nyq);
        set(handles.statustext,'String','Computing currv');
        pause(.01);
        handles.currv=cfunct(handles.theta,handles.zeta,handles.data.currvmnc,handles.data.xm_nyq,xn_nyq);
        set(handles.statustext,'String','Computing dR/du');
        pause(.01);
        handles.drdu=sfunct(handles.theta,handles.zeta,handles.data.rumns,handles.data.xm,xn);
        set(handles.statustext,'String','Computing dR/dv');
        pause(.01);
        handles.drdv=sfunct(handles.theta,handles.zeta,handles.data.rvmns,handles.data.xm,xn);
        set(handles.statustext,'String','Computing dZ/du');
        pause(.01);
        handles.dzdu=cfunct(handles.theta,handles.zeta,handles.data.zumnc,handles.data.xm,xn);
        set(handles.statustext,'String','Computing dZ/dv');
        pause(.01);
        handles.dzdv=cfunct(handles.theta,handles.zeta,handles.data.zvmnc,handles.data.xm,xn);
        if isfield(handles,'prot')
            handles = rmfield(handles,'prot');
            handles = rmfield(handles,'U_rot');
            handles = rmfield(handles,'prpr');
            handles = rmfield(handles,'vphi');
        end
        if isfield(handles.data,'protmnc')
            pause(.01);
            handles.prot=cfunct(handles.theta,handles.zeta,handles.data.protmnc,handles.data.xm,xn);
            pause(.01);
            handles.U_rot=cfunct(handles.theta,handles.zeta,handles.data.protrsqmnc,handles.data.xm,xn);
            pause(.01);
            handles.prpr=cfunct(handles.theta,handles.zeta,handles.data.prprmnc,handles.data.xm,xn);
        end
        % Handle asymetric vars
        if handles.data.iasym==1
            if handles.diag_plt==1
                handles.rbs=permute(handles.data.rbs, [3 1 2]);
                handles.zbc=permute(handles.data.zbc, [3 1 2]);
            end
            set(handles.statustext,'String','Asymetric VMEC detected!');
            pause(.01);
            set(handles.statustext,'String','Adding rbs to R');
            pause(.01);
            handles.r=handles.r+sfunct(handles.theta,handles.zeta,handles.data.rmns,handles.data.xm,xn);
            handles.drdu=handles.drdu+cfunct(handles.theta,handles.zeta,handles.data.rumnc,handles.data.xm,xn);
            handles.drdv=handles.drdv+cfunct(handles.theta,handles.zeta,handles.data.rvmnc,handles.data.xm,xn);
            set(handles.statustext,'String','Adding zbc to Z');
            pause(.01);
            handles.z=handles.z+cfunct(handles.theta,handles.zeta,handles.data.zmnc,handles.data.xm,xn);
            handles.dzdu=handles.dzdu+sfunct(handles.theta,handles.zeta,handles.data.zumns,handles.data.xm,xn);
            handles.dzdv=handles.dzdv+sfunct(handles.theta,handles.zeta,handles.data.zvmns,handles.data.xm,xn);
            handles.lam=handles.lam+cfunct(handles.theta,handles.zeta,handles.data.lmns,handles.data.xm,xn);
            % Handle new asymetric VMEC values
            if isfield(handles.data,'bmns')
                set(handles.statustext,'String','Adding bs to b');
                pause(.01);
                handles.b=handles.b+sfunct(handles.theta,handles.zeta,handles.data.bmns,handles.data.xm_nyq,xn_nyq);
            end
            if isfield(handles.data,'gmns')
                set(handles.statustext,'String','Adding gs to g');
                pause(.01);
                handles.g=handles.g+sfunct(handles.theta,handles.zeta,handles.data.gmns,handles.data.xm_nyq,xn_nyq);
            end
            if isfield(handles.data,'bsubsmnc')
                set(handles.statustext,'String','Adding b_sc to b_s');
                pause(.01);
                handles.b_s=handles.b_s+cfunct(handles.theta,handles.zeta,handles.data.bsubsmnc,handles.data.xm_nyq,xn_nyq);
            end
            if isfield(handles.data,'bsubumns')
                set(handles.statustext,'String','Adding b_us to b_u');
                pause(.01);
                handles.b_u=handles.b_u+sfunct(handles.theta,handles.zeta,handles.data.bsubumns,handles.data.xm_nyq,xn_nyq);
            end
            if isfield(handles.data,'bsubvmns')
                set(handles.statustext,'String','Adding b_vs to b_v');
                pause(.01);
                handles.b_v=handles.b_v+sfunct(handles.theta,handles.zeta,handles.data.bsubvmns,handles.data.xm_nyq,xn_nyq);
            end
            if isfield(handles.data,'currvmns')
                set(handles.statustext,'String','Adding currvs to currv');
                pause(.01);
                handles.currv=handles.currv+sfunct(handles.theta,handles.zeta,handles.data.currvmns,handles.data.xm_nyq,xn_nyq);
                %handles.b=handles.b+sfunct(handles.theta,handles.zeta,handles.data.bmns,handles.data.xm_nyq,xn);
            end
            if isfield(handles.data,'protmns')
                pause(.01);
                handles.prot=handles.prot+sfunct(handles.theta,handles.zeta,handles.data.protmns,handles.data.xm,xn);
                pause(.01);
                handles.U_rot=handles.U_rot+sfunct(handles.theta,handles.zeta,handles.data.protrsqmns,handles.data.xm,xn);
                pause(.01);
                handles.prpr=handles.prpr+sfunct(handles.theta,handles.zeta,handles.data.prprmns,handles.data.xm,xn);
            end
        end
        % Fix jcurv as they are multiplied by gsqrt
        handles.currv=handles.currv./handles.g;
        % Calculate cylindrical components of B
        handles.br=handles.bu.*handles.drdu+handles.bv.*handles.drdv;
        handles.bphi=handles.r.*handles.bv;
        handles.bz=handles.bu.*handles.dzdu+handles.bv.*handles.dzdv;
        % Caculate poloidal components of B
        handles.br_pol=handles.bu.*handles.drdu;
        handles.bz_pol=handles.bu.*handles.dzdu;
        % Add in the toroidal current
        handles.jphi=handles.r.*handles.currv;
        % Calculate the toroidal velocity
        if isfield(handles.data,'omega')
            handles.data.vphi = handles.r.*0.0;
            %VA = handles.b(1,1,1)/sqrt(pi.*4E-7.*handles.prot(1,1,1).^(1/handles.data.gamma));
            %mass = pi*4E-7*handles.data.presf(1)*abs(handles.data.vp(1))^handles.data.gamma;
            %cs = sqrt(2*handles.data.presf(1)/mass);
            %cs = sqrt(abs(handles.data.vp(1)).^(1-handles.data.gamma)./(pi*4E-7));
            cs = sqrt(handles.data.tpotb(1).*1.60217733E-19./1.6726231E-27); % Assume [eV] and proton mass
            for j=1:handles.data.ns
               handles.vphi(j,:,:) = cs.*sqrt(handles.data.machsq).*handles.data.omega(j).*handles.r(j,:,:)./handles.data.Rmajor;
            end
        end
        if isfield(handles.data,'currumnc')
            set(handles.statustext,'String','Computing curru');
            pause(.01);
            handles.curru=cfunct(handles.theta,handles.zeta,handles.data.currumnc,handles.data.xm_nyq,xn_nyq);
            if handles.data.iasym==1
                if isfield(handles.data,'currumns')
                    handles.curru=handles.curru+sfunct(handles.theta,handles.zeta,handles.data.currumns,handles.data.xm_nyq,xn_nyq);
                end
            end
            handles.curru=handles.curru./handles.g;
            handles.jr=handles.curru.*handles.drdu+handles.currv.*handles.drdv;
            handles.jz=handles.curru.*handles.dzdu+handles.currv.*handles.dzdv;
            handles.jpara=handles.jr.*handles.br+handles.jphi.*handles.bphi+handles.jz.*handles.bz;
            handles.jpara=handles.jpara./handles.b;
            if handles.diag_plt==1
                handles.rbc=permute(handles.data.rbc, [3 1 2]);
                handles.zbs=permute(handles.data.zbs, [3 1 2]);
                set(handles.statustext,'String','Computing FORCE NORM');
                pause(.01);
                handles.jx=0.0.*handles.jr;
                handles.jy=0.0.*handles.jr;
                handles.jz=0.0.*handles.jr;
                for j=1:handles.ntor+1
                    handles.jx(:,:,j)=handles.jr(:,:,j).*cos(handles.zeta(j))-handles.jphi(:,:,j).*sin(handles.zeta(j));
                    handles.jy(:,:,j)=handles.jr(:,:,j).*sin(handles.zeta(j))+handles.jphi(:,:,j).*cos(handles.zeta(j));
                    handles.bx(:,:,j)=handles.br(:,:,j).*cos(handles.zeta(j))-handles.bphi(:,:,j).*sin(handles.zeta(j));
                    handles.by(:,:,j)=handles.br(:,:,j).*sin(handles.zeta(j))+handles.bphi(:,:,j).*cos(handles.zeta(j));
                end
                set(handles.statustext,'String','Computing JxB');
                pause(.01);
                handles.jxbx=handles.jy.*handles.bz-handles.jz.*handles.by;
                handles.jxby=handles.jz.*handles.bx-handles.jx.*handles.bz;
                handles.jxbz=handles.jx.*handles.by-handles.jy.*handles.bx;
                set(handles.statustext,'String','Computing grad(p)');
                pause(.01);
                handles.p=0.0*handles.jxbx;
                handles.dpds=0.0*handles.jxbx;
                for i=1:handles.data.ns
                    handles.p(i,:,:)=handles.data.presf(i);
                end
                handles.dpds(1,:,:)=0.0;
                handles.dpds(handles.data.ns,:,:)=handles.p(handles.data.ns,:,:)-handles.p(handles.data.ns-1,:,:);
                for i=2:handles.data.ns-1
                    handles.dpds(i,:,:)=0.5*(handles.p(i+1,:,:)-handles.p(i-1,:,:));
                end
                handles.dpds=handles.dpds.*handles.data.ns;
                % Calculate d/ds (s=i/ns)
                drds=0.0*handles.jxbx;
                dzds=0.0*handles.jxbx;
                for i=2:handles.data.ns-1
                    drds(i,:,:)=0.5*(handles.r(i+1,:,:)-handles.r(i-1,:,:));
                    dzds(i,:,:)=0.5*(handles.z(i+1,:,:)-handles.z(i-1,:,:));
                end
                drds(1,:,:)=handles.r(2,:,:)-handles.r(1,:,:);
                drds(handles.data.ns,:,:)=handles.r(handles.data.ns,:,:)-handles.r(handles.data.ns-1,:,:);
                dzds(1,:,:)=handles.z(2,:,:)-handles.z(1,:,:);
                dzds(handles.data.ns,:,:)=handles.z(handles.data.ns,:,:)-handles.z(handles.data.ns-1,:,:);
                drds=drds.*handles.data.ns;
                dzds=dzds.*handles.data.ns;
                set(handles.statustext,'String','Computing Fn');
                pause(.01);
                % Now compute Fp
                fpr=handles.dpds.*drds;
                fpz=handles.dpds.*dzds;
                for j=1:handles.ntor+1
                    fpx(:,:,j)=fpr(:,:,j).*cos(handles.zeta(j));
                    fpy(:,:,j)=fpr(:,:,j).*sin(handles.zeta(j));
                end
                handles.fnx=handles.jxbx-fpx;
                handles.fny=handles.jxby-fpy;
                handles.fnz=handles.jxbz-fpz;
                handles.fn=sqrt(handles.fnx.*handles.fnx+handles.fny.*handles.fny+handles.fnz.*handles.fnz);
            end
        else % Remove the jr and jz handles if they exist
            if isfield(handles,'jr'), rmfield(handles,'jr'); end
            if isfield(handles,'jz'), rmfield(handles,'jz'); end
        end
    case 'SPEC'% Reset theta to 1D
        handles.xlabel='Normalize Toroidal Flux (s)';
        handles.theta=0:2*pi/(handles.mpol-1):2*pi;
        % Transform quantities
        set(handles.statustext,'String','Computing r');
        pause(.01);
        handles.r=cfunct(handles.theta,handles.zeta,handles.data.rmnc,handles.data.xm,handles.data.xn);
        if isfield(handles.data,'rmns')
            handles.r=handles.r+sfunct(handles.theta,handles.zeta,handles.data.rmns,handles.data.xm,handles.data.xn);
        end
        set(handles.statustext,'String','Computing z');
        pause(.01);
        handles.z=sfunct(handles.theta,handles.zeta,handles.data.zmns,handles.data.xm,handles.data.xn);
        if isfield(handles.data,'zmnc')
            handles.z=handles.z+cfunct(handles.theta,handles.zeta,handles.data.zmnc,handles.data.xm,handles.data.xn);
        end
        set(handles.statustext,'String','Computing B-Field');
        pause(.01);
        %%%% (mu-nv)
        for i = 1:handles.data.ns
            rumns(:,i)=-handles.data.rmnc(:,i).*handles.data.xm';
            rumnc(:,i)=handles.data.rmns(:,i).*handles.data.xm';
            zumns(:,i)=-handles.data.zmnc(:,i).*handles.data.xm';
            zumnc(:,i)=handles.data.zmns(:,i).*handles.data.xm';
            rvmns(:,i)=handles.data.rmnc(:,i).*handles.data.xn';
            rvmnc(:,i)=-handles.data.rmns(:,i).*handles.data.xn';
            zvmns(:,i)=handles.data.zmnc(:,i).*handles.data.xn';
            zvmnc(:,i)=-handles.data.zmns(:,i).*handles.data.xn';
        end
        ru=cfunct(handles.theta,handles.zeta,rumnc,handles.data.xm,handles.data.xn);
        ru=ru+sfunct(handles.theta,handles.zeta,rumns,handles.data.xm,handles.data.xn);
        zu=cfunct(handles.theta,handles.zeta,zumnc,handles.data.xm,handles.data.xn);
        zu=zu+sfunct(handles.theta,handles.zeta,zumns,handles.data.xm,handles.data.xn);
        rs=cfunct(handles.theta,handles.zeta,handles.data.rmncp,handles.data.xm,handles.data.xn);
        rs=rs+sfunct(handles.theta,handles.zeta,handles.data.rmnsp,handles.data.xm,handles.data.xn);
        zs=cfunct(handles.theta,handles.zeta,handles.data.zmncp,handles.data.xm,handles.data.xn);
        zs=zs+sfunct(handles.theta,handles.zeta,handles.data.zmnsp,handles.data.xm,handles.data.xn);
        handles.g = handles.r.*(ru.*zs-rs.*zu);
        handles.g(1,:,:) = 2.0.*handles.g(2,:,:)-handles.g(3,:,:);
        %handles.g(1,:,:) = 0.50*handles.r(1,:,:); % Asympotitic Jacobian value on axis
        rv=cfunct(handles.theta,handles.zeta,rvmnc,handles.data.xm,handles.data.xn);
        rv=rv+sfunct(handles.theta,handles.zeta,rvmns,handles.data.xm,handles.data.xn);
        zv=cfunct(handles.theta,handles.zeta,zvmnc,handles.data.xm,handles.data.xn);
        zv=zv+sfunct(handles.theta,handles.zeta,zvmns,handles.data.xm,handles.data.xn);
        for mn = 1:handles.data.mnmax
            bsupsmnc(mn,:) =  ( handles.data.xm(mn).*handles.data.Asubvmns(mn,:) + handles.data.xn(mn).*handles.data.Asubumns(mn,:));
            bsupsmns(mn,:) = -( handles.data.xm(mn).*handles.data.Asubvmnc(mn,:) + handles.data.xn(mn).*handles.data.Asubumnc(mn,:));
            bsupumnc(mn,:) = -handles.data.Asubvmncp(mn,:);
            bsupumns(mn,:) = -handles.data.Asubvmnsp(mn,:);
            bsupvmnc(mn,:) =  handles.data.Asubumncp(mn,:);
            bsupvmns(mn,:) =  handles.data.Asubumnsp(mn,:);
        end
        handles.bs=cfunct(handles.theta,handles.zeta,bsupsmnc,handles.data.xm,handles.data.xn);
        handles.bs=handles.bs+sfunct(handles.theta,handles.zeta,bsupsmns,handles.data.xm,handles.data.xn);
        handles.bs=handles.bs./handles.g;
        handles.bu=cfunct(handles.theta,handles.zeta,bsupumnc,handles.data.xm,handles.data.xn);
        handles.bu=handles.bu+sfunct(handles.theta,handles.zeta,bsupumns,handles.data.xm,handles.data.xn);
        handles.bu=handles.bu./handles.g;
        handles.bv=cfunct(handles.theta,handles.zeta,bsupvmnc,handles.data.xm,handles.data.xn);
        handles.bv=handles.bv+sfunct(handles.theta,handles.zeta,bsupvmns,handles.data.xm,handles.data.xn);
        handles.bv=handles.bv./handles.g;
        handles.br=handles.bu.*ru+handles.bv.*rv;
        handles.bphi=handles.r.*handles.bv;
        handles.bz=handles.bu.*zu+handles.bv.*zv;
        handles.ju=[]; handles.jv=[]; handles.jr=[]; handles.jphi=[]; handles.jz=[];
        for i=1:handles.data.ns
            handles.ju(i,:,:)   = handles.bu(i,:,:).*handles.data.muf(i)./(4*pi*1E-7);
            handles.jv(i,:,:)   = handles.bv(i,:,:).*handles.data.muf(i)./(4*pi*1E-7);
            handles.jr(i,:,:)   = handles.br(i,:,:).*handles.data.muf(i)./(4*pi*1E-7);
            handles.jphi(i,:,:) = handles.bphi(i,:,:).*handles.data.muf(i)./(4*pi*1E-7);
            handles.jz(i,:,:)   = handles.bz(i,:,:).*handles.data.muf(i)./(4*pi*1E-7);
        end
        handles.data.jcuru = sum(sum(handles.ju,3),2)./(4*pi*pi);
        handles.data.jcurv = sum(sum(handles.jv,3),2)./(4*pi*pi);
    case 'siesta'
        handles.xlabel='Normalized minor radius (r/a)';
        handles.r=cfunct(handles.theta,handles.zeta,handles.data.rmnc,handles.data.xm,handles.data.xn);
        handles.z=sfunct(handles.theta,handles.zeta,handles.data.zmns,handles.data.xm,handles.data.xn);
        handles.p=cfunct(handles.theta,handles.zeta,handles.data.pmnc,handles.data.xm,handles.data.xn);
        handles.bs=sfunct(handles.theta,handles.zeta,handles.data.bsupsmns,handles.data.xm,handles.data.xn);
        handles.bu=cfunct(handles.theta,handles.zeta,handles.data.bsupumnc,handles.data.xm,handles.data.xn);
        handles.bv=cfunct(handles.theta,handles.zeta,handles.data.bsupvmnc,handles.data.xm,handles.data.xn);
        handles.b_s=sfunct(handles.theta,handles.zeta,handles.data.bsubsmns,handles.data.xm,handles.data.xn);
        handles.b_u=cfunct(handles.theta,handles.zeta,handles.data.bsubumnc,handles.data.xm,handles.data.xn);
        handles.b_v=cfunct(handles.theta,handles.zeta,handles.data.bsubvmnc,handles.data.xm,handles.data.xn);
        handles.currs=sfunct(handles.theta,handles.zeta,handles.data.currsmns,handles.data.xm,handles.data.xn);
        handles.curru=cfunct(handles.theta,handles.zeta,handles.data.currumnc,handles.data.xm,handles.data.xn);
        handles.currv=cfunct(handles.theta,handles.zeta,handles.data.currvmnc,handles.data.xm,handles.data.xn);
        handles.bu=cfunct(handles.theta,handles.zeta,handles.data.bsupumnc,handles.data.xm,handles.data.xn);
        handles.bv=cfunct(handles.theta,handles.zeta,handles.data.bsupvmnc,handles.data.xm,handles.data.xn);
        handles.drds=cfunct(handles.theta,handles.zeta,handles.data.rsmnc,handles.data.xm,handles.data.xn);
        handles.drdu=sfunct(handles.theta,handles.zeta,handles.data.rumns,handles.data.xm,handles.data.xn);
        handles.drdv=sfunct(handles.theta,handles.zeta,handles.data.rvmns,handles.data.xm,handles.data.xn);
        handles.dzds=sfunct(handles.theta,handles.zeta,handles.data.zsmns,handles.data.xm,handles.data.xn);
        handles.dzdu=cfunct(handles.theta,handles.zeta,handles.data.zumnc,handles.data.xm,handles.data.xn);
        handles.dzdv=cfunct(handles.theta,handles.zeta,handles.data.zvmnc,handles.data.xm,handles.data.xn);
        % Calculate jacbian
        handles.g=handles.r.*(handles.drdu.*handles.dzds-handles.drds.*handles.dzdu);
        % Fix current
        handles.currs = handles.currs./handles.g;
        handles.cvrru = handles.curru./handles.g;
        handles.currv = handles.currv./handles.g;
        % Calc |B|
        handles.b = sqrt(handles.bs.*handles.bs+handles.bu.*handles.bu+handles.bv.*handles.bv);
        % Calculate cylindrical components of B
        handles.br=handles.bs.*handles.drds+handles.bu.*handles.drdu+handles.bv.*handles.drdv;
        handles.bphi=handles.r.*handles.bv;
        handles.bz=handles.bs.*handles.dzds+handles.bu.*handles.dzdu+handles.bv.*handles.dzdv;
        % Calculate cylindrical components of J
        handles.jr=handles.currs.*handles.drds+handles.curru.*handles.drdu+handles.currv.*handles.drdv;
        handles.jphi=handles.r.*handles.currv;
        handles.jz=handles.currs.*handles.dzds+handles.curru.*handles.dzdu+handles.currv.*handles.dzdv;
    case 'nout'% NSTAB DATA
        handles.xlabel='Normalize Toroidal Flux (s)';
        handles.theta=handles.data.theta;
        handles.zeta=handles.data.zeta;
        %handles.r = cfunct(handles.theta,handles.zeta,handles.data.rhomn,handles.data.xm,handles.data.xn);
        %handles.z = sfunct(handles.theta,handles.zeta,handles.data.rhomn,handles.data.xm,handles.data.xn);
        %handles.r = handles.r + cfunct(handles.theta,handles.zeta,handles.data.r0mnc,handles.data.xm,handles.data.xn);
        %handles.z = handles.z + sfunct(handles.theta,handles.zeta,handles.data.z0mns,handles.data.xm,handles.data.xn);
        %for i=2:handles.data.ns-1
        %    handles.r(i,:,:) = (1.0-handles.r(i,:,:)).*handles.r(1,:,:) + handles.r(i,:,:) .* handles.r(handles.data.ns,:,:);
        %    handles.z(i,:,:) = (1.0-handles.r(i,:,:)).*handles.z(1,:,:) + handles.r(i,:,:) .* handles.z(handles.data.ns,:,:);
        %end
        %handles.r = handles.r*handles.data.nstabscale+handles.data.r0_vmec;
        %handles.z = handles.z*handles.data.nstabscale;
        handles.r=handles.data.rreal;
        handles.z=handles.data.zreal;
    case 'pies_out'
        handles.xlabel='Normalize Toroidal Flux (s)';
        % Note for the derivatives terms we assume (mu+nv) but PIES (nv-mu)
        % This is handled in the read_pies_netcdf routine but not for the
        % derivative terms so the du drivative terms we must multiply by a
        % negative number.  We pass negative theta to the functions to
        % handle the PIES kernel.
        set(handles.statustext,'String','Computing r');
        pause(.01);
        handles.r=cfunct(-handles.theta,handles.zeta,handles.data.rmnc,handles.data.xm,handles.data.xn);
        set(handles.statustext,'String','Computing z');
        pause(.01);
        handles.z=sfunct(-handles.theta,handles.zeta,handles.data.zmns,handles.data.xm,handles.data.xn);
        set(handles.statustext,'String','Computing Brho');
        pause(.01);
        handles.bs=sfunct(-handles.theta,handles.zeta,handles.data.bsmnc,handles.data.xm,handles.data.xn);
        set(handles.statustext,'String','Computing Btheta');
        pause(.01);
        handles.bu=cfunct(-handles.theta,handles.zeta,handles.data.bumnc,handles.data.xm,handles.data.xn);
        set(handles.statustext,'String','Computing Bphi');
        pause(.01);
        handles.bv=cfunct(-handles.theta,handles.zeta,handles.data.bvmnc,handles.data.xm,handles.data.xn);
        if isfield(handles.data,'jlc')
            set(handles.statustext,'String','Computing J_||');
            pause(.01);
            handles.jpara=cfunct(-handles.theta,handles.zeta,handles.data.jlmnc,handles.data.xm,handles.data.xn);
        end
        if isfield(handles.data,'press')
            set(handles.statustext,'String','Computing Pressure');
            pause(.01);
            handles.press=cfunct(-handles.theta,handles.zeta,handles.data.pressmnc,handles.data.xm,handles.data.xn);
        end
        set(handles.statustext,'String','Computing dR/du');
        pause(.01);
        handles.drdu=-sfunct(-handles.theta,handles.zeta,handles.data.rumns,handles.data.xm,handles.data.xn);
        set(handles.statustext,'String','Computing dR/dv');
        pause(.01);
        handles.drdv=sfunct(-handles.theta,handles.zeta,handles.data.rvmns,handles.data.xm,handles.data.xn);
        set(handles.statustext,'String','Computing dZ/du');
        pause(.01);
        handles.dzdu=-cfunct(-handles.theta,handles.zeta,handles.data.zumnc,handles.data.xm,handles.data.xn);
        set(handles.statustext,'String','Computing dZ/dv');
        pause(.01);
        handles.dzdv=cfunct(-handles.theta,handles.zeta,handles.data.zvmnc,handles.data.xm,handles.data.xn);
        % Calculate d/ds (s=i/ns)
        for i=2:handles.data.ns-1
            drds(i,:,:)=0.5*(handles.r(i+1,:,:)-handles.r(i-1,:,:));
            dzds(i,:,:)=0.5*(handles.z(i+1,:,:)-handles.z(i-1,:,:));
        end
        drds(1,:,:)=handles.r(2,:,:)-handles.r(1,:,:);
        drds(handles.data.ns,:,:)=handles.r(handles.data.ns,:,:)-handles.r(handles.data.ns-1,:,:);
        dzds(1,:,:)=handles.z(2,:,:)-handles.z(1,:,:);
        dzds(handles.data.ns,:,:)=handles.z(handles.data.ns,:,:)-handles.z(handles.data.ns-1,:,:);
        drds=drds./handles.data.ns;
        dzds=dzds./handles.data.ns;
        % Compute 
        handles.bphi=handles.r.*handles.bv;
        handles.br=handles.bs.*drds+handles.bu.*handles.drdu+handles.bv.*handles.drdv;
        handles.bz=handles.bs.*dzds+handles.bu.*handles.dzdu+handles.bv.*handles.dzdv;
        % Mod B
        handles.b=sqrt(handles.br.*handles.br+handles.bphi.*handles.bphi+handles.bz.*handles.bz);
        % B-NORM
        snr=handles.dzdu.*handles.r;
        snphi=handles.drdu.*handles.dzdv-handles.drdv.*handles.dzdu;
        %for i=1:handles.ntor
        %    snx(:,:,i)=snr(:,:,i).*cos(handles.zeta(i))-snphi(:,:,i).*sin(handles.zeta(i));
        %    sny(:,:,i)=snr(:,:,i).*sin(handles.zeta(i))+snphi(:,:,i).*cos(handles.zeta(i));
        %end
        snz=-handles.drdu.*handles.r;
        sn=sqrt(snr.*snr+snphi.*snphi+snz.*snz);
        handles.snr=snr./sn;
        handles.snphi=snphi./sn;
        handles.snz=snz./sn;
        handles.bnorm=handles.br.*handles.snr+handles.bphi.*handles.snphi+...
            handles.bz.*handles.snz;
        if handles.diag_plt==1
            handles.rbc=permute(handles.data.rbc, [3 1 2]);
            handles.zbs=permute(handles.data.zbs, [3 1 2]);
        end
    case 'boozer'   % For Boozer Coordinates transform
        handles.xlabel='Normalize Toroidal Flux (s)';
        % Reduce boozer quantities in radial dimension
        handles.data.rbc=handles.data.rbc(:,:,handles.data.idx==1);
        handles.data.zbs=handles.data.zbs(:,:,handles.data.idx==1);
        handles.data.ps=handles.data.ps(:,:,handles.data.idx==1);
        handles.data.gc=handles.data.gc(:,:,handles.data.idx==1);
        handles.data.bc=handles.data.bc(:,:,handles.data.idx==1);
        if handles.data.lasym
            handles.data.rbs=handles.data.rbs(:,:,handles.data.idx==1);
            handles.data.zbc=handles.data.zbc(:,:,handles.data.idx==1);
            handles.data.pc=handles.data.pc(:,:,handles.data.idx==1);
            handles.data.gs=handles.data.gs(:,:,handles.data.idx==1);
            handles.data.bs=handles.data.bs(:,:,handles.data.idx==1);
        end
        %handles.data.ns=sum(handles.data.idx);
        % Reset theta and zeta to 1D
        handles.theta=0:2*pi/(handles.mpol-1):2*pi;
        handles.zeta=0:2*pi/(handles.ntor-1):2*pi;
        pause(.01);
        % Calculate zeta transform
        set(handles.statustext,'String','Computing ptran');
        pause(.01);
        handles.ptran=sfunct(handles.theta,handles.zeta,handles.data.pmns,handles.data.xm,handles.data.xn);
        if handles.data.lasym
            handles.ptran=handles.ptran+...
                cfunct(handles.theta,handles.zeta,handles.data.pmnc,handles.data.xm,handles.data.xn);
        end
        % Transform quantities
        set(handles.statustext,'String','Computing r');
        pause(.01);
        handles.r=cfunct(handles.theta,handles.zeta,handles.data.rmnc,handles.data.xm,handles.data.xn);
        set(handles.statustext,'String','Computing z');
        pause(.01);
        handles.z=sfunct(handles.theta,handles.zeta,handles.data.zmns,handles.data.xm,handles.data.xn);
        set(handles.statustext,'String','Computing |B|');
        pause(.01);
        handles.b=cfunct(handles.theta,handles.zeta,handles.data.bmnc,handles.data.xm,handles.data.xn);
        set(handles.statustext,'String','Computing g');
        pause(.01);
        handles.g=cfunct(handles.theta,handles.zeta,handles.data.gmnc,handles.data.xm,handles.data.xn);
        if handles.data.lasym
            set(handles.statustext,'String','Computing r-sin');
            pause(.01);
            handles.r=handles.r+...
                sfunct(handles.theta,handles.zeta,handles.data.rmns,handles.data.xm,handles.data.xn);
            set(handles.statustext,'String','Computing z-cos');
            pause(.01);
            handles.z=handles.z+...
                cfunct(handles.theta,handles.zeta,handles.data.zmnc,handles.data.xm,handles.data.xn);
            set(handles.statustext,'String','Computing |B|-sin');
            pause(.01);
            handles.b=handles.b+...
                sfunct(handles.theta,handles.zeta,handles.data.bmns,handles.data.xm,handles.data.xn);
            set(handles.statustext,'String','Computing g-sin');
            pause(.01);
            handles.g=hangles.g+...
                sfunct(handles.theta,handles.zeta,handles.data.gmns,handles.data.xm,handles.data.xn);
        end
    case 'FIELDLINES'
        handles.r = handles.data.raxis;
        handles.z = handles.data.zaxis;
        nlines = handles.data.nlines;
        for j=1:handles.data.npoinc
            temp = handles.data.R_lines(1:nlines,j:handles.data.npoinc:handles.data.nsteps);
            temp_R(1:nlines,1:size(temp,2),j) = temp;
            temp = handles.data.Z_lines(1:nlines,j:handles.data.npoinc:handles.data.nsteps);
            temp_Z(1:nlines,1:size(temp,2),j) = temp;
        end
        handles.data.R_lines = permute(temp_R,[1 3 2]);
        handles.data.Z_lines = permute(temp_Z,[1 3 2]);
    case 'BEAM3D'
        handles.data.ntor = handles.data.nphi;
    otherwise
        set(handles.statustext,'String','Datatype:ERROR');
        text(0.35,0.5,['ERROR: Unknown Datatype: ' handles.data.datatype],...
            'Color','red');
        return
end
% Done With transformations
set(handles.statustext,'String','Ready');
handles.xmax=max(max(max(handles.r)));
handles.xmin=min(min(min(handles.r)));
handles.ymax=handles.xmax;
handles.ymin=handles.xmin;
handles.zmax=max(max(max(handles.z)));
handles.zmin=min(min(min(handles.z)));
guidata(hObject,handles);


% --- Executes when selected object is changed in cutplane.
function cutplane_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in cutplane 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag')
    case 'rcut'
        handles.cuttype='r1d';
    case 'thetacut'
        handles.cuttype='theta1d';
    case 'zetacut'
        handles.cuttype='zeta1d';
    case 'fluxcut'
        handles.cuttype='r2';
    case 'polcut'
        handles.cuttype='theta2';
    case 'rzcut'
        handles.cuttype='zeta2';
    case 'threedcut'
        handles.cuttype='3D';
end
guidata(hObject,handles);
guiupdate(handles);


% --- Executes on slider movement.
function theslide_Callback(hObject, eventdata, handles)
% hObject    handle to theslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.thetaval=int32(get(hObject,'Value'));
guidata(hObject,handles);
update_plots(handles);

% --- Executes during object creation, after setting all properties.
function theslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function rslide_Callback(hObject, eventdata, handles)
% hObject    handle to rslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.rval=int32(get(hObject,'Value'));
guidata(hObject,handles);
update_plots(handles);

% --- Executes during object creation, after setting all properties.
function rslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function outfile_Callback(hObject, eventdata, handles)
% hObject    handle to outfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outfile as text
%        str2double(get(hObject,'String')) returns contents of outfile as a double


% --- Executes during object creation, after setting all properties.
function outfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=get(handles.outfile,'String');
c = strfind(filename,'.');
c = max(c);
ddex=c+1;
%[S c]=textscan(filename,'%[^.]');
%ddex=c+2;
nlen=size(filename,2);
lvrml=0;
if ddex < nlen
    filename(ddex:nlen)=lower(filename(ddex:nlen));
    switch filename(ddex:nlen)
        case {'ai','bmp','emf','eps','fig','jpg','m','pbm','pcx','pdf',...
                'ps','pgm','png','ppm','tif'}
            ext_temp = filename(ddex:nlen);
            file_temp = filename(1:ddex-2);
        case {'wrl'}
            lvrml=1;
        otherwise
            filename=strcat(filename(1:c),'.fig');
            set(handles.statustext,'String','Saving as matlab fig!','ForegroundColor','red');
            pause(1.0);
    end
else
    filename=strcat(filename,'.fig');
end
h=figure('Visible','on','Renderer','opengl');
update_plots(handles);
drawnow;
if lvrml
    vrml(h,filename);
else
    switch ext_temp
        case{'ps'}
            saveas(h,file_temp,'psc');
            
        case{'eps'}
            saveas(h,file_temp,'epsc');
                
        otherwise
            saveas(h,filename);
    end
end
set(handles.statustext,'String','Ready','ForegroundColor','black');
drawnow;
close(h);
