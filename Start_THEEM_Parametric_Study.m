function varargout = Start_THEEM_Parametric_Study(varargin)
% Start_THEEM_Parametric_Study MATLAB code for Start_THEEM_Parametric_Study.fig
%      Start_THEEM_Parametric_Study, by itself, creates a new Start_THEEM_Parametric_Study or raises the existing
%      singleton*.
%
%      H = Start_THEEM_Parametric_Study returns the handle to a new Start_THEEM_Parametric_Study or the handle to
%      the existing singleton*.
%
%      Start_THEEM_Parametric_Study('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Start_THEEM_Parametric_Study.M with the given input arguments.
%      Start_THEEM_Parametric_Study('Property','Value',...) creates a new Start_THEEM_Parametric_Study or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Start_THEEM_Parametric_Study_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Start_THEEM_Parametric_Study_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Start_THEEM_Parametric_Study

% Last Modified by GUIDE v2.5 18-Feb-2016 15:21:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Start_THEEM_Parametric_Study_OpeningFcn, ...
                   'gui_OutputFcn',  @Start_THEEM_Parametric_Study_OutputFcn, ...
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


% --- Executes just before Start_THEEM_Parametric_Study is made visible.
function Start_THEEM_Parametric_Study_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Start_THEEM_Parametric_Study (see VARARGIN)

% Choose default command line output for Start_THEEM_Parametric_Study
handles.output = hObject;
% handles.Cancel_program='Cancel';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Start_THEEM_Parametric_Study wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Start_THEEM_Parametric_Study_OutputFcn(hObject, eventdata, handles) 
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
set(handles.figure1, 'HandleVisibility', 'off'); %Turns off access to GUI. Exempts GUI from next line.
close all; %Close all open figures except for GUI
set(handles.figure1, 'HandleVisibility', 'on'); %Turns on access to GUI again.
clc;evalin('base','clear'); %Clears command window & base workspace
[FileName,PathName] = uigetfile('*.mat'); %User selects existing input file
uiresume
close(handles.figure1); 
if FileName~=0 %A .mat file needs to be selected
   vars = whos('-file',[PathName,FileName]); %Stores variables from .mat file without loading into workspace
   mat_size=numel(vars); %Number of variables in mat file
   if mat_size==27 %Check to see if mat file has the correct number of inputs to run THEEM
       load([PathName,FileName],'THEEM_model') %Specify which THEEM model is being used, i.e. 2-D vs. 3-D
       compare=strcmp([PathName,FileName],[pwd,'\Optimization Program\Parametric Study\THEEM_Input_Parametric.mat']); %Compares string of selected file to string of default input file
       if compare==0 %If file is existing .mat file in folder, program does not copy it.  It just uses it.
       copyfile([PathName,FileName],[pwd,'\Optimization Program\Parametric Study\THEEM_Input_Parametric.mat'])
       end
       CTGH_Parametric_Study
   else %If incorrect input file, will give error and restart program
       errordlg({'This is not a Parametric Study input file.','Please select another file.'},'Input File Error')
       uiwait(gcf);
       Start_THEEM_Parametric_Study
   end
else
   Start_THEEM_Parametric_Study %Restarts program if no .mat file is selected
end



% --- Executes on button press in New_button.
function New_button_Callback(hObject, eventdata, handles)
% hObject    handle to New_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global THEEM_model
THEEM_model='Parametric Study';
varargout=CTGH_Input_GUI({THEEM_model});
cancel='Cancel';
if strcmp(varargout,cancel)==0    
    guidata(hObject, handles);
    uiresume
    close all; %Close all open figures except for GUI
    clc;evalin('base','clear'); %Clears command window & base workspace
    load('THEEM_Input_Parametric.mat','THEEM_model'); 
    CTGH_Parametric_Study;
%     evalin('base','load(''THEEM_Output_Parametric.mat'')');
end




% --- Executes on button press in Cancel_button.
function Cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.Cancel_program='Cancel';
guidata(hObject, handles);
uiresume
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
end
