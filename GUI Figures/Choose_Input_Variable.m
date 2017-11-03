function varargout = Choose_Input_Variable(varargin)
% CHOOSE_INPUT_VARIABLE MATLAB code for Choose_Input_Variable.fig
%      CHOOSE_INPUT_VARIABLE, by itself, creates a new CHOOSE_INPUT_VARIABLE or raises the existing
%      singleton*.
%
%      H = CHOOSE_INPUT_VARIABLE returns the handle to a new CHOOSE_INPUT_VARIABLE or the handle to
%      the existing singleton*.
%
%      CHOOSE_INPUT_VARIABLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHOOSE_INPUT_VARIABLE.M with the given input arguments.
%
%      CHOOSE_INPUT_VARIABLE('Property','Value',...) creates a new CHOOSE_INPUT_VARIABLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Choose_Input_Variable_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Choose_Input_Variable_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Choose_Input_Variable

% Last Modified by GUIDE v2.5 26-Oct-2017 16:22:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Choose_Input_Variable_OpeningFcn, ...
                   'gui_OutputFcn',  @Choose_Input_Variable_OutputFcn, ...
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


% --- Executes just before Choose_Input_Variable is made visible.
function Choose_Input_Variable_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Choose_Input_Variable (see VARARGIN)

% Choose default command line output for Choose_Input_Variable
handles.output = hObject;
handles.Cancel_program=char('True');
% inputs=fieldnames(load('THEEM_Input_Parametric.mat'));
inputs=varargin{1};
set(handles.Variable_selection,'string',inputs)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Choose_Input_Variable wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Choose_Input_Variable_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.Cancel_program;
if isfield(handles, 'input_variable')~=0 %Checks to see if user selected a variable
    varargout{2} = handles.input_variable;
    varargout{3} = handles.input_min_value;
    varargout{4} = handles.input_max_value;
    varargout{5} = handles.input_step_value;
else
    varargout{2} = [];
    varargout{3} = [];
    varargout{4} = [];
    varargout{5} = [];
end
 delete(hObject);




% --- Executes on button press in Cancel_Button.
function Cancel_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume
guidata(hObject, handles);



% --- Executes on selection change in Variable_selection.
function Variable_selection_Callback(hObject, eventdata, handles)
% hObject    handle to Variable_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Variable_selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Variable_selection


% --- Executes during object creation, after setting all properties.
function Variable_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Variable_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Run_Button.
function Run_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Run_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.input_min_value=str2num(char(get(handles.Min_Value,'String')));
handles.input_max_value=str2num(char(get(handles.Max_Value,'String')));
handles.input_step_value=str2num(char(get(handles.Step_Value,'String')));
input_index=get(handles.Variable_selection,'Value');
input_variable_list=get(handles.Variable_selection,'String');
handles.input_variable=char(input_variable_list(input_index));
handles.Cancel_program=char('False');
uiresume
guidata(hObject, handles); 



function Min_Value_Callback(hObject, eventdata, handles)
% hObject    handle to Min_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Min_Value as text
%        str2double(get(hObject,'String')) returns contents of Min_Value as a double


% --- Executes during object creation, after setting all properties.
function Min_Value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Min_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Max_Value_Callback(hObject, eventdata, handles)
% hObject    handle to Max_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') %returns contents of Max_Value as text
%        str2double(get(hObject,'String')) returns contents of Max_Value as a double


% --- Executes during object creation, after setting all properties.
function Max_Value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Max_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Step_Value_Callback(hObject, eventdata, handles)
% hObject    handle to Step_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Step_Value as text
%        str2double(get(hObject,'String')) returns contents of Step_Value as a double


% --- Executes during object creation, after setting all properties.
function Step_Value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Step_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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
