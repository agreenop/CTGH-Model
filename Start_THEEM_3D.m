function varargout = Start_THEEM_3D(varargin)
% START_THEEM MATLAB code for Start_THEEM_3D.fig
%      START_THEEM, by itself, creates a new START_THEEM or raises the existing
%      singleton*.
%
%      H = START_THEEM returns the handle to a new START_THEEM or the handle to
%      the existing singleton*.
%
%      START_THEEM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in START_THEEM.M with the given input arguments.
%
%      START_THEEM('Property','Value',...) creates a new START_THEEM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Start_THEEM_3D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Start_THEEM_3D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Start_THEEM_3D

% Last Modified by GUIDE v2.5 18-Feb-2016 15:21:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Start_THEEM_3D_OpeningFcn, ...
                   'gui_OutputFcn',  @Start_THEEM_3D_OutputFcn, ...
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


% --- Executes just before Start_THEEM_3D is made visible.
function Start_THEEM_3D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Start_THEEM_3D (see VARARGIN)

% Choose default command line output for Start_THEEM_3D
handles.output = hObject;
% handles.Cancel_program='Cancel';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Start_THEEM_3D wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Start_THEEM_3D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.Cancel_program;
delete(hObject);


% --- Executes on button press in Existing_button.
function Existing_button_Callback(hObject, eventdata, handles)
% hObject    handle to Existing_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.Cancel_program='Run';
guidata(hObject, handles);
close(handles.figure1); 
load('THEEM_3D_Input.mat');
if isequal(model_selection,'Test Bundle 1')
    clc;clear;
    run('Mockup1_2D.m')
    load('THEEM_Input.mat');
    load('THEEM_Output.mat');
else
    clc;clear;
    run('CTGH_3D.m')
    load('THEEM_3D_Input.mat');
    load('THEEM_3D_Output.mat');
end

% --- Executes on button press in New_button.
function New_button_Callback(hObject, eventdata, handles)
% hObject    handle to New_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CTGH_Input_GUI;
load('THEEM_3D_Input.mat');
% handles.Cancel_program=program;
guidata(hObject, handles);
close(handles.figure1);
if isequal(model_selection,'Test Bundle 1')
    clc;clear;
    run('Mockup1_2D.m')
    load('THEEM_Input.mat');
    load('THEEM_Output.mat');
else
    clc;clear;
    run('CTGH_3D.m')
    load('THEEM_3D_Input.mat');
    load('THEEM_3D_Output.mat');
end




% --- Executes on button press in Cancel_button.
function Cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.Cancel_program='Cancel';
guidata(hObject, handles);
close(handles.figure1); 
return;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
% delete(hObject);
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, call UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
    load('THEEM_Output.mat');
end
