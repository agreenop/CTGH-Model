function varargout = CTGH_Input_GUI(varargin)
% CTGH_Input_GUI MATLAB code for CTGH_Input_GUI.fig
%      CTGH_Input_GUI, by itself, creates a new CTGH_Input_GUI or raises the existing
%      singleton*.
%
%      H = CTGH_Input_GUI returns the handle to a new CTGH_Input_GUI or the handle to
%      the existing singleton*.
%
%      CTGH_Input_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CTGH_Input_GUI.M with the given input arguments.
%
%      CTGH_Input_GUI('Property','Value',...) creates a new CTGH_Input_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CTGH_Input_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CTGH_Input_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CTGH_Input_GUI

% Last Modified by GUIDE v2.5 21-Apr-2016 15:58:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CTGH_Input_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CTGH_Input_GUI_OutputFcn, ...
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


% --- Executes just before CTGH_Input_GUI is made visible.
function CTGH_Input_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CTGH_Input_GUI (see VARARGIN)

% Choose default command line output for CTGH_Input_GUI
handles.output = hObject;
handles.Cancel_program='Cancel';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CTGH_Input_GUI wait for user response (see UIRESUME)
uiwait(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = CTGH_Input_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.Cancel_program;
 delete(hObject);


% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gas_index=get(handles.Gas_Type,'Value');
gas_list=get(handles.Gas_Type,'String');
gas=char(gas_list(gas_index));
liqud_index=get(handles.Liquid_Type,'Value');
liquid_list=get(handles.Liquid_Type,'String');
liquid=char(liquid_list(liqud_index));
T_g_in=str2num(get(handles.Inlet_Gas_Temp,'String'));
T_l_in=str2num(get(handles.Inlet_Liquid_Temp,'String'));
P_g_in=str2num(get(handles.Inlet_Gas_Press,'String'));
P_l_in=str2num(get(handles.Inlet_Liquid_Press,'String'));
m_g=str2num(get(handles.Gas_Mass_Flow,'String'));
m_l=str2num(get(handles.Liquid_Mass_Flow,'String'));
material_index=get(handles.Tube_Material,'Value');
material_list=cellstr(get(handles.Tube_Material,'String'));
tube_material=char(material_list(material_index));
units_list=get(handles.Diameter_Units,'String');
diameter_index=get(handles.Diameter_Units,'Value');
diameter_units=char(units_list(diameter_index));
D_out=str2num(get(handles.Tube_Diameter,'String'));
switch diameter_units
    case 'in'
        D_out=D_out*0.0254;
    case 'ft'
         D_out=D_out*0.0254/12;
    case 'mm'
         D_out=D_out*1000;
    case 'cm'
         D_out=D_out*100;
    case 'm'
         D_out=D_out;
end
thick_index=get(handles.Thickness_Units,'Value');
thick_units=char(units_list(thick_index));
t=str2num(get(handles.Tube_Thickness,'String'));
switch thick_units
    case 'in'
         t=t*0.0254;
    case 'ft'
         t=t*0.0254/12;
    case 'mm'
         t=t*1000;
    case 'cm'
         t=t*100;
    case 'm'
        t=t;
end
SL=str2num(get(handles.Long_Pitch_Ratio,'String'));
ST=str2num(get(handles.Transverse_Pitch_Ratio,'String'));
entry=str2num(get(handles.Liquid_Inlets,'String'));
model_index=get(handles.Model_Selection,'Value');
model_list=get(handles.Model_Selection,'String');
model_selection=char(model_list(model_index));
handles.Cancel_program='Run';
save('THEEM_3D_Input.mat','gas','liquid','T_g_in','T_l_in','P_g_in','P_l_in','m_g','m_l','tube_material','D_out','t','SL','ST','entry','model_selection')
guidata(hObject, handles);
close(handles.figure1); 

% --- Executes on selection change in Gas_Type.
function Gas_Type_Callback(hObject, eventdata, handles)
% hObject    handle to Gas_Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Gas_Type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Gas_Type


% --- Executes during object creation, after setting all properties.
function Gas_Type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gas_Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Liquid_Type.
function Liquid_Type_Callback(hObject, eventdata, handles)
% hObject    handle to Liquid_Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Liquid_Type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Liquid_Type


% --- Executes during object creation, after setting all properties.
function Liquid_Type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Liquid_Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Inlet_Gas_Temp_Callback(hObject, eventdata, handles)
% hObject    handle to Inlet_Gas_Temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Inlet_Gas_Temp as text
%        str2double(get(hObject,'String')) returns contents of Inlet_Gas_Temp as a double


% --- Executes during object creation, after setting all properties.
function Inlet_Gas_Temp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Inlet_Gas_Temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Inlet_Gas_Press_Callback(hObject, eventdata, handles)
% hObject    handle to Inlet_Gas_Press (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Inlet_Gas_Press as text
%        str2double(get(hObject,'String')) returns contents of Inlet_Gas_Press as a double


% --- Executes during object creation, after setting all properties.
function Inlet_Gas_Press_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Inlet_Gas_Press (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Inlet_Liquid_Temp_Callback(hObject, eventdata, handles)
% hObject    handle to Liquid_Mass_Flow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Liquid_Mass_Flow as text
%        str2double(get(hObject,'String')) returns contents of Liquid_Mass_Flow as a double


% --- Executes during object creation, after setting all properties.
function Inlet_Liquid_Temp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Liquid_Mass_Flow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Inlet_Liquid_Press_Callback(hObject, eventdata, handles)
% hObject    handle to Inlet_Liquid_Temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Inlet_Liquid_Temp as text
%        str2double(get(hObject,'String')) returns contents of Inlet_Liquid_Temp as a double


% --- Executes during object creation, after setting all properties.
function Inlet_Liquid_Press_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Inlet_Liquid_Temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Liquid_Mass_Flow_Callback(hObject, eventdata, handles)
% hObject    handle to Inlet_Liquid_Press (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Inlet_Liquid_Press as text
%        str2double(get(hObject,'String')) returns contents of Inlet_Liquid_Press as a double


% --- Executes during object creation, after setting all properties.
function Liquid_Mass_Flow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Inlet_Liquid_Press (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gas_Mass_Flow_Callback(hObject, eventdata, handles)
% hObject    handle to Liquid_Mass_Flow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Liquid_Mass_Flow as text
%        str2double(get(hObject,'String')) returns contents of Liquid_Mass_Flow as a double


% --- Executes during object creation, after setting all properties.
function Gas_Mass_Flow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Liquid_Mass_Flow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Default.
function Default_Callback(hObject, eventdata, handles)
% hObject    handle to Default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model_index=get(handles.Model_Selection,'Value');
model_list=get(handles.Model_Selection,'String');
Model_Selection=char(model_list(model_index));
switch Model_Selection
    case 'Mk1 CTGH'
set(handles.Gas_Type,'Value',1);
set(handles.Liquid_Type,'Value',1);
set(handles.Tube_Material,'Value',1);
set(handles.Diameter_Units,'Value',1);
set(handles.Thickness_Units,'Value',1);
set(handles.Inlet_Gas_Temp,'String',418.6);
set(handles.Inlet_Gas_Press,'String',18.76);
set(handles.Gas_Mass_Flow,'String',418.5);
set(handles.Inlet_Liquid_Temp,'String',700);
set(handles.Inlet_Liquid_Press,'String',2);
set(handles.Liquid_Mass_Flow,'String',480.2);
set(handles.Tube_Diameter,'String',0.25);
set(handles.Tube_Thickness,'String',0.035);
set(handles.Long_Pitch_Ratio,'String',1.45);
set(handles.Transverse_Pitch_Ratio,'String',1.256);
set(handles.Liquid_Inlets,'String',4);
    case 'Test Bundle 1'
set(handles.Gas_Type,'Value',1);
set(handles.Liquid_Type,'Value',2);
set(handles.Tube_Material,'Value',1);
set(handles.Diameter_Units,'Value',1);
set(handles.Thickness_Units,'Value',1);
set(handles.Inlet_Gas_Temp,'String',23.7);
set(handles.Inlet_Gas_Press,'String',1);
set(handles.Gas_Mass_Flow,'String',.34);
set(handles.Inlet_Liquid_Temp,'String',49.7);
set(handles.Inlet_Liquid_Press,'String',2);
set(handles.Liquid_Mass_Flow,'String',0.2);
set(handles.Tube_Diameter,'String',0.25);
set(handles.Tube_Thickness,'String',0.002);
set(handles.Long_Pitch_Ratio,'String',1.37);
set(handles.Transverse_Pitch_Ratio,'String',1.67);
set(handles.Liquid_Inlets,'String',2);
end

% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Cancel_program='Cancel';
guidata(hObject, handles);
close(handles.figure1); 
return;



function Liquid_Inlets_Callback(hObject, eventdata, handles)
% hObject    handle to Liquid_Inlets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Liquid_Inlets as text
%        str2double(get(hObject,'String')) returns contents of Liquid_Inlets as a double


% --- Executes during object creation, after setting all properties.
function Liquid_Inlets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Liquid_Inlets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Transverse_Pitch_Ratio_Callback(hObject, eventdata, handles)
% hObject    handle to Transverse_Pitch_Ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Transverse_Pitch_Ratio as text
%        str2double(get(hObject,'String')) returns contents of Transverse_Pitch_Ratio as a double


% --- Executes during object creation, after setting all properties.
function Transverse_Pitch_Ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Transverse_Pitch_Ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Long_Pitch_Ratio_Callback(hObject, eventdata, handles)
% hObject    handle to Long_Pitch_Ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Long_Pitch_Ratio as text
%        str2double(get(hObject,'String')) returns contents of Long_Pitch_Ratio as a double


% --- Executes during object creation, after setting all properties.
function Long_Pitch_Ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Long_Pitch_Ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tube_Thickness_Callback(hObject, eventdata, handles)
% hObject    handle to Tube_Thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tube_Thickness as text
%        str2double(get(hObject,'String')) returns contents of Tube_Thickness as a double


% --- Executes during object creation, after setting all properties.
function Tube_Thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tube_Thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tube_Diameter_Callback(hObject, eventdata, handles)
% hObject    handle to Tube_Diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tube_Diameter as text
%        str2double(get(hObject,'String')) returns contents of Tube_Diameter as a double


% --- Executes during object creation, after setting all properties.
function Tube_Diameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tube_Diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Tube_Material.
function Tube_Material_Callback(hObject, eventdata, handles)
% hObject    handle to Tube_Material (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Tube_Material contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Tube_Material


% --- Executes during object creation, after setting all properties.
function Tube_Material_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tube_Material (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Diameter_Units.
function Diameter_Units_Callback(hObject, eventdata, handles)
% hObject    handle to Diameter_Units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Diameter_Units contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Diameter_Units


% --- Executes during object creation, after setting all properties.
function Diameter_Units_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Diameter_Units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Thickness_Units.
function Thickness_Units_Callback(hObject, eventdata, handles)
% hObject    handle to Thickness_Units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Thickness_Units contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Thickness_Units


% --- Executes during object creation, after setting all properties.
function Thickness_Units_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Thickness_Units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in Model_Selection.
function Model_Selection_Callback(hObject, eventdata, handles)
% hObject    handle to Model_Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model_index=get(handles.Model_Selection,'Value');
model_list=get(handles.Model_Selection,'String');
Model_Selection=char(model_list(model_index));
switch Model_Selection
    case 'Mk1 CTGH'
set(handles.Gas_Type,'Value',1);
set(handles.Liquid_Type,'Value',1);
set(handles.Tube_Material,'Value',1);
set(handles.Diameter_Units,'Value',1);
set(handles.Thickness_Units,'Value',1);
set(handles.Inlet_Gas_Temp,'String',418.6);
set(handles.Inlet_Gas_Press,'String',18.76);
set(handles.Gas_Mass_Flow,'String',418.5);
set(handles.Inlet_Liquid_Temp,'String',700);
set(handles.Inlet_Liquid_Press,'String',2);
set(handles.Liquid_Mass_Flow,'String',480.2);
set(handles.Tube_Diameter,'String',0.25);
set(handles.Tube_Thickness,'String',0.035);
set(handles.Long_Pitch_Ratio,'String',1.45);
set(handles.Transverse_Pitch_Ratio,'String',1.256);
set(handles.Liquid_Inlets,'String',4);
    case 'Test Bundle 1'
set(handles.Gas_Type,'Value',1);
set(handles.Liquid_Type,'Value',2);
set(handles.Tube_Material,'Value',1);
set(handles.Diameter_Units,'Value',1);
set(handles.Thickness_Units,'Value',1);
set(handles.Inlet_Gas_Temp,'String',23.7);
set(handles.Inlet_Gas_Press,'String',1);
set(handles.Gas_Mass_Flow,'String',.34);
set(handles.Inlet_Liquid_Temp,'String',49.7);
set(handles.Inlet_Liquid_Press,'String',2);
set(handles.Liquid_Mass_Flow,'String',0.2);
set(handles.Tube_Diameter,'String',0.25);
set(handles.Tube_Thickness,'String',0.02);
set(handles.Long_Pitch_Ratio,'String',1.37);
set(handles.Transverse_Pitch_Ratio,'String',1.67);
set(handles.Liquid_Inlets,'String',2);
end

% Hints: contents = cellstr(get(hObject,'String')) returns Model_Selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Model_Selection


% --- Executes during object creation, after setting all properties.
function Model_Selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Model_Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
