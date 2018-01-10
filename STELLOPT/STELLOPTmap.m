function varargout = STELLOPTmap(varargin)
% STELLOPTMAP MATLAB code for STELLOPTmap.fig
%      STELLOPTMAP, by itself, creates a new STELLOPTMAP or raises the existing
%      singleton*.
%
%      H = STELLOPTMAP returns the handle to a new STELLOPTMAP or the handle to
%      the existing singleton*.
%
%      STELLOPTMAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STELLOPTMAP.M with the given input arguments.
%
%      STELLOPTMAP('Property','Value',...) creates a new STELLOPTMAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before STELLOPTmap_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to STELLOPTmap_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help STELLOPTmap

% Last Modified by GUIDE v2.5 05-Feb-2013 13:48:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @STELLOPTmap_OpeningFcn, ...
                   'gui_OutputFcn',  @STELLOPTmap_OutputFcn, ...
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


% --- Executes just before STELLOPTmap is made visible.
function STELLOPTmap_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to STELLOPTmap (see VARARGIN)

% Choose default command line output for STELLOPTmap
handles.output = hObject;
handles.data = read_stellopt('map.dat');
set(handles.xaxis_menu,'String',handles.data.var_name);
set(handles.xaxis_menu,'Value',1.0);
set(handles.yaxis_menu,'String',handles.data.var_name);
set(handles.yaxis_menu,'String',2.0);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes STELLOPTmap wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = STELLOPTmap_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in yaxis_menu.
function yaxis_menu_Callback(hObject, eventdata, handles)
% hObject    handle to yaxis_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns yaxis_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from yaxis_menu
dex_x = get(handles.xaxis_menu,'Value');
dex_y = get(handles.yaxis_menu,'Value');
index_x = handles.data.x(:,dex_x)==handles.data.x(1,dex_x);



% --- Executes during object creation, after setting all properties.
function yaxis_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yaxis_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in xaxis_menu.
function xaxis_menu_Callback(hObject, eventdata, handles)
% hObject    handle to xaxis_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xaxis_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xaxis_menu


% --- Executes during object creation, after setting all properties.
function xaxis_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xaxis_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
