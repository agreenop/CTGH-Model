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

% Last Modified by GUIDE v2.5 21-Oct-2016 17:13:01

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
THEEM_model=varargin{1}{1};
handles.THEEM_model=THEEM_model;
% Update handles structure
guidata(hObject, handles);
% get(handles.CTGH_Input_GUI,'UserData');
switch THEEM_model
    case '0D'
        set(handles.tube_holders,'Enable','off')
    case '2D'
        set(handles.Outlet_Liquid_Temp,'Enable','off')
        set(handles.Outlet_Gas_Temp,'Enable','off')
    case '3D'
        set(handles.Outlet_Liquid_Temp,'Enable','off')
        set(handles.Outlet_Gas_Temp,'Enable','off')        
    case 'Optimization'
         set(handles.Tube_Diameter,'Enable','off')
         set(handles.Diameter_Units,'Enable','off')
         set(handles.Tube_Thickness,'Enable','off')
         set(handles.Thickness_Units,'Enable','off')
         set(handles.Long_Pitch_Ratio,'Enable','off')
         set(handles.Transverse_Pitch_Ratio,'Enable','off')
         set(handles.Liquid_Inlets,'Enable','off')
         set(handles.Loops,'Enable','off')
         set(handles.tube_layer,'Enable','off')
         set(handles.layer_num,'Enable','off')
         set(handles.bundles,'Enable','off')
         set(handles.spacers,'Enable','off')
         set(handles.spacer_width,'Enable','off')
         set(handles.tube_holders,'Enable','off')
         set(handles.Bundle_Radius,'Enable','off')
         set(handles.Gap_Units,'Enable','off')
         set(handles.Bundle_Radius_Units,'Enable','off')
         set(handles.tube_slope,'Enable','off')
         set(handles.heat_rod,'Enable','off')
        end
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
function Run_Callback(hObject, eventdata, handles, varargin)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
THEEM_model=handles.THEEM_model;
gas_index=get(handles.Gas_Type,'Value');
gas_list=get(handles.Gas_Type,'String');
gas=char(gas_list(gas_index));
liqud_index=get(handles.Liquid_Type,'Value');
liquid_list=get(handles.Liquid_Type,'String');
liquid=char(liquid_list(liqud_index));
T_g_in=str2num(get(handles.Inlet_Gas_Temp,'String'));
T_g_out=str2num(get(handles.Outlet_Gas_Temp,'String'));
T_l_in=str2num(get(handles.Inlet_Liquid_Temp,'String'));
T_l_out=str2num(get(handles.Outlet_Liquid_Temp,'String'));
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
         D_out=D_out/1000;
    case 'cm'
         D_out=D_out/100;
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
         t=t/1000;
    case 'cm'
         t=t/100;
    case 'm'
        t=t;
end
SL=str2num(get(handles.Long_Pitch_Ratio,'String'));
ST=str2num(get(handles.Transverse_Pitch_Ratio,'String'));
entry=str2num(get(handles.Liquid_Inlets,'String'));
model_index=get(handles.Model_Selection,'Value');
model_list=get(handles.Model_Selection,'String');
model_selection=char(model_list(model_index));
loops=str2num(get(handles.Loops,'String'));
tube_layer=str2num(get(handles.tube_layer,'String'));
layer_num=str2num(get(handles.layer_num,'String'));
bundles=str2num(get(handles.bundles,'String'));
spacers=str2num(get(handles.spacers,'String'));
gap_width_index=get(handles.Gap_Units,'Value');
gap_units_list=get(handles.Gap_Units,'String');
gap_units=char(gap_units_list(gap_width_index));
spacer_width=str2num(get(handles.spacer_width,'String'));
switch gap_units
    case 'in'
        spacer_width=spacer_width*0.0254;
    case 'ft'
         spacer_width=spacer_width*0.0254/12;
    case 'mm'
         spacer_width=spacer_width/1000;
    case 'cm'
         spacer_width=spacer_width/100;
    case 'm'
         spacer_width=spacer_width;
end
tube_holders=str2num(get(handles.tube_holders,'String'));
bundle_radius_index=get(handles.Bundle_Radius_Units,'Value');
bundle_radius_units_list=get(handles.Bundle_Radius_Units,'String');
radius_units=char(bundle_radius_units_list(bundle_radius_index));
R_ci=str2num(get(handles.Bundle_Radius,'String'));
switch radius_units
    case 'in'
        R_ci=R_ci*0.0254;
    case 'ft'
         R_ci=R_ci*0.0254/12;
    case 'mm'
         R_ci=R_ci/1000;
    case 'cm'
         R_ci=R_ci/100;
    case 'm'
         R_ci=R_ci;
end
tube_slope=str2num(get(handles.tube_slope,'String'));
heat_rod=str2num(get(handles.heat_rod,'String'));
handles.Cancel_program='Run';
switch THEEM_model
    case '0D'
        save('0-D Model/THEEM_Input_0D.mat','gas','liquid','T_g_in','T_l_in',...
            'T_g_out','T_l_out','P_g_in','P_l_in','m_g','m_l','tube_material',...
            'D_out','t','SL','ST','entry','THEEM_model','loops','tube_layer',...
            'layer_num','bundles','spacers','spacer_width','R_ci','tube_slope','heat_rod')
    case '2D'
        save('2-D Model/THEEM_Input_2D.mat','gas','liquid','T_g_in','T_l_in','P_g_in','P_l_in',...
            'm_g','m_l','tube_material','D_out','t','SL','ST','entry','THEEM_model','loops',...
            'tube_layer','layer_num','bundles','spacers','spacer_width','tube_holders',...
            'R_ci','tube_slope','heat_rod')
    case '3D'
        save('3-D Model/THEEM_Input_3D.mat','gas','liquid','T_g_in','T_l_in','P_g_in','P_l_in',...
            'm_g','m_l','tube_material','D_out','t','SL','ST','entry','THEEM_model','loops',...
            'tube_layer','layer_num','bundles','spacers','spacer_width','tube_holders',...
            'R_ci','tube_slope','heat_rod')        
    case 'Optimization'
        save('Optimization Program/THEEM_Input_Optimization.mat','gas','liquid',...
            'T_g_in','T_l_in','T_g_out','T_l_out','P_g_in','P_l_in','m_g',...
            'm_l','tube_material')
        end    
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


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



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
        set(handles.Outlet_Gas_Temp,'String',670);
        set(handles.Inlet_Gas_Press,'String',18.76);
        set(handles.Gas_Mass_Flow,'String',418.5);
        set(handles.Inlet_Liquid_Temp,'String',700);
        set(handles.Outlet_Liquid_Temp,'String',600);
        set(handles.Inlet_Liquid_Press,'String',3.5);
        set(handles.Liquid_Mass_Flow,'String',480.2);
        set(handles.Tube_Diameter,'String',0.25);
        set(handles.Tube_Thickness,'String',0.035);
        set(handles.Long_Pitch_Ratio,'String',1.256);
        set(handles.Transverse_Pitch_Ratio,'String',1.45);
        set(handles.Liquid_Inlets,'String',4);
        set(handles.Loops,'String',3);
        set(handles.tube_layer,'String',3);
        set(handles.layer_num,'String',40);
        set(handles.bundles,'String',36);
        set(handles.spacers,'String',2);
        set(handles.spacer_width,'String',0.038);
        set(handles.Gap_Units,'Value',1);
        set(handles.tube_holders,'String',12);
        set(handles.Bundle_Radius,'String',0.662);
        set(handles.Bundle_Radius_Units,'Value',1);
        set(handles.tube_slope,'String',0.003);
        set(handles.heat_rod,'String',1/2);
    case 'Test Bundle 1'
        set(handles.Gas_Type,'Value',1);
        set(handles.Liquid_Type,'Value',2);
        set(handles.Tube_Material,'Value',1);
        set(handles.Diameter_Units,'Value',1);
        set(handles.Thickness_Units,'Value',1);
        set(handles.Inlet_Gas_Temp,'String',23.7);
        set(handles.Outlet_Gas_Temp,'String',43);
        set(handles.Inlet_Gas_Press,'String',1);
        set(handles.Gas_Mass_Flow,'String',.34);
        set(handles.Inlet_Liquid_Temp,'String',49.7);
        set(handles.Outlet_Liquid_Temp,'String',41);
        set(handles.Inlet_Liquid_Press,'String',2);
        set(handles.Liquid_Mass_Flow,'String',0.2);
        set(handles.Tube_Diameter,'String',0.25);
        set(handles.Tube_Thickness,'String',0.02);
        set(handles.Long_Pitch_Ratio,'String',1.37);
        set(handles.Transverse_Pitch_Ratio,'String',1.67);
        set(handles.Liquid_Inlets,'String',2);
        set(handles.Loops,'String',4);
        set(handles.tube_layer,'String',1);
        set(handles.layer_num,'String',20);
        set(handles.bundles,'String',1);
        set(handles.spacers,'String',0);
        set(handles.spacer_width,'String',0);
        set(handles.Gap_Units,'Value',1);
        set(handles.tube_holders,'String',6);
        set(handles.Bundle_Radius,'String',8.75);
        set(handles.Bundle_Radius_Units,'Value',4);
        set(handles.tube_slope,'String',0);
        set(handles.heat_rod,'String',0);
    case 'CASET'
        set(handles.Gas_Type,'Value',1);
        set(handles.Liquid_Type,'Value',2);
        set(handles.Tube_Material,'Value',1);
        set(handles.Diameter_Units,'Value',1);
        set(handles.Thickness_Units,'Value',1);
        set(handles.Inlet_Gas_Temp,'String',25);
        set(handles.Outlet_Gas_Temp,'String',51);
        set(handles.Inlet_Gas_Press,'String',1);
        set(handles.Gas_Mass_Flow,'String',0.649);
        set(handles.Inlet_Liquid_Temp,'String',80);
        set(handles.Outlet_Liquid_Temp,'String',40);
        set(handles.Inlet_Liquid_Press,'String',2);
        set(handles.Liquid_Mass_Flow,'String',0.1);
        set(handles.Tube_Diameter,'String',0.25);
        set(handles.Tube_Thickness,'String',0.02);
        set(handles.Long_Pitch_Ratio,'String',1.256);
        set(handles.Transverse_Pitch_Ratio,'String',1.50);
        set(handles.Liquid_Inlets,'String',2);
        set(handles.Loops,'String',3);
        set(handles.tube_layer,'String',2);
        set(handles.layer_num,'String',10);
        set(handles.bundles,'String',1);
        set(handles.spacers,'String',1);
        set(handles.spacer_width,'String',0.038);
        set(handles.Gap_Units,'Value',1);
        set(handles.tube_holders,'String',6);
        set(handles.Bundle_Radius,'String',25);
        set(handles.Bundle_Radius_Units,'Value',2);
        set(handles.tube_slope,'String',0);
        set(handles.heat_rod,'String',0);
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

function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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
model_selection=char(model_list(model_index));
switch model_selection
    case 'Mk1 CTGH'
        set(handles.Gas_Type,'Value',1);
        set(handles.Liquid_Type,'Value',1);
        set(handles.Tube_Material,'Value',1);
        set(handles.Diameter_Units,'Value',1);
        set(handles.Thickness_Units,'Value',1);
        set(handles.Inlet_Gas_Temp,'String',418.6);
        set(handles.Outlet_Gas_Temp,'String',670);
        set(handles.Inlet_Gas_Press,'String',18.76);
        set(handles.Gas_Mass_Flow,'String',418.5);
        set(handles.Inlet_Liquid_Temp,'String',700);
        set(handles.Outlet_Liquid_Temp,'String',600);
        set(handles.Inlet_Liquid_Press,'String',3.5);
        set(handles.Liquid_Mass_Flow,'String',480.2);
        set(handles.Tube_Diameter,'String',0.25);
        set(handles.Tube_Thickness,'String',0.035);
        set(handles.Long_Pitch_Ratio,'String',1.256);
        set(handles.Transverse_Pitch_Ratio,'String',1.45);
        set(handles.Liquid_Inlets,'String',4);
        set(handles.Loops,'String',3);
        set(handles.tube_layer,'String',3);
        set(handles.layer_num,'String',40);
        set(handles.bundles,'String',36);
        set(handles.spacers,'String',2);
        set(handles.spacer_width,'String',0.038);
        set(handles.Gap_Units,'Value',1);
        set(handles.tube_holders,'String',12);
        set(handles.Bundle_Radius,'String',0.662);
        set(handles.Bundle_Radius_Units,'Value',1);
        set(handles.tube_slope,'String',0.003);
        set(handles.heat_rod,'String',1/2);
    case 'Test Bundle 1'
        set(handles.Gas_Type,'Value',1);
        set(handles.Liquid_Type,'Value',2);
        set(handles.Tube_Material,'Value',1);
        set(handles.Diameter_Units,'Value',1);
        set(handles.Thickness_Units,'Value',1);
        set(handles.Inlet_Gas_Temp,'String',23.7);
        set(handles.Outlet_Gas_Temp,'String',43);
        set(handles.Inlet_Gas_Press,'String',1);
        set(handles.Gas_Mass_Flow,'String',.34);
        set(handles.Inlet_Liquid_Temp,'String',49.7);
        set(handles.Outlet_Liquid_Temp,'String',41);
        set(handles.Inlet_Liquid_Press,'String',2);
        set(handles.Liquid_Mass_Flow,'String',0.2);
        set(handles.Tube_Diameter,'String',0.25);
        set(handles.Tube_Thickness,'String',0.02);
        set(handles.Long_Pitch_Ratio,'String',1.37);
        set(handles.Transverse_Pitch_Ratio,'String',1.67);
        set(handles.Liquid_Inlets,'String',2);
        set(handles.Loops,'String',4);
        set(handles.tube_layer,'String',1);
        set(handles.layer_num,'String',20);
        set(handles.bundles,'String',1);
        set(handles.spacers,'String',0);
        set(handles.spacer_width,'String',0);
        set(handles.Gap_Units,'Value',1);
        set(handles.tube_holders,'String',6);
        set(handles.Bundle_Radius,'String',8.75);
        set(handles.Bundle_Radius_Units,'Value',4);
        set(handles.tube_slope,'String',0);
        set(handles.heat_rod,'String',0);
    case 'CASET'
        set(handles.Gas_Type,'Value',1);
        set(handles.Liquid_Type,'Value',2);
        set(handles.Tube_Material,'Value',1);
        set(handles.Diameter_Units,'Value',1);
        set(handles.Thickness_Units,'Value',1);
        set(handles.Inlet_Gas_Temp,'String',25);
        set(handles.Outlet_Gas_Temp,'String',51);
        set(handles.Inlet_Gas_Press,'String',1);
        set(handles.Gas_Mass_Flow,'String',0.649);
        set(handles.Inlet_Liquid_Temp,'String',80);
        set(handles.Outlet_Liquid_Temp,'String',40);
        set(handles.Inlet_Liquid_Press,'String',2);
        set(handles.Liquid_Mass_Flow,'String',0.1);
        set(handles.Tube_Diameter,'String',0.25);
        set(handles.Tube_Thickness,'String',0.02);
        set(handles.Long_Pitch_Ratio,'String',1.256);
        set(handles.Transverse_Pitch_Ratio,'String',1.50);
        set(handles.Liquid_Inlets,'String',2);
        set(handles.Loops,'String',3);
        set(handles.tube_layer,'String',2);
        set(handles.layer_num,'String',10);
        set(handles.bundles,'String',1);
        set(handles.spacers,'String',1);
        set(handles.spacer_width,'String',0.038);
        set(handles.Gap_Units,'Value',1);
        set(handles.tube_holders,'String',6);
        set(handles.Bundle_Radius,'String',25);
        set(handles.Bundle_Radius_Units,'Value',2);
        set(handles.tube_slope,'String',0);
        set(handles.heat_rod,'String',0);
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



function Loops_Callback(hObject, eventdata, handles)
% hObject    handle to Loops (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Loops as text
%        str2double(get(hObject,'String')) returns contents of Loops as a double


% --- Executes during object creation, after setting all properties.
function Loops_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Loops (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tube_layer_Callback(hObject, eventdata, handles)
% hObject    handle to tube_layer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tube_layer as text
%        str2double(get(hObject,'String')) returns contents of tube_layer as a double


% --- Executes during object creation, after setting all properties.
function tube_layer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tube_layer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function layer_num_Callback(hObject, eventdata, handles)
% hObject    handle to layer_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of layer_num as text
%        str2double(get(hObject,'String')) returns contents of layer_num as a double


% --- Executes during object creation, after setting all properties.
function layer_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to layer_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bundles_Callback(hObject, eventdata, handles)
% hObject    handle to bundles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bundles as text
%        str2double(get(hObject,'String')) returns contents of bundles as a double


% --- Executes during object creation, after setting all properties.
function bundles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bundles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spacers_Callback(hObject, eventdata, handles)
% hObject    handle to spacers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spacers as text
%        str2double(get(hObject,'String')) returns contents of spacers as a double


% --- Executes during object creation, after setting all properties.
function spacers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spacers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spacer_width_Callback(hObject, eventdata, handles)
% hObject    handle to spacer_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spacer_width as text
%        str2double(get(hObject,'String')) returns contents of spacer_width as a double


% --- Executes during object creation, after setting all properties.
function spacer_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spacer_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tube_holders_Callback(hObject, eventdata, handles)
% hObject    handle to tube_holders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tube_holders as text
%        str2double(get(hObject,'String')) returns contents of tube_holders as a double


% --- Executes during object creation, after setting all properties.
function tube_holders_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tube_holders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bundle_Radius_Callback(hObject, eventdata, handles)
% hObject    handle to Bundle_Radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bundle_Radius as text
%        str2double(get(hObject,'String')) returns contents of Bundle_Radius as a double


% --- Executes during object creation, after setting all properties.
function Bundle_Radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bundle_Radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Gap_Units.
function Gap_Units_Callback(hObject, eventdata, handles)
% hObject    handle to Gap_Units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Gap_Units contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Gap_Units


% --- Executes during object creation, after setting all properties.
function Gap_Units_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gap_Units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Bundle_Radius_Units.
function Bundle_Radius_Units_Callback(hObject, eventdata, handles)
% hObject    handle to Bundle_Radius_Units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Bundle_Radius_Units contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Bundle_Radius_Units


% --- Executes during object creation, after setting all properties.
function Bundle_Radius_Units_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bundle_Radius_Units (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tube_slope_Callback(hObject, eventdata, handles)
% hObject    handle to tube_slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tube_slope as text
%        str2double(get(hObject,'String')) returns contents of tube_slope as a double


% --- Executes during object creation, after setting all properties.
function tube_slope_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tube_slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function heat_rod_Callback(hObject, eventdata, handles)
% hObject    handle to heat_rod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of heat_rod as text
%        str2double(get(hObject,'String')) returns contents of heat_rod as a double


% --- Executes during object creation, after setting all properties.
function heat_rod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heat_rod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Outlet_Gas_Temp_Callback(hObject, eventdata, handles)
% hObject    handle to Outlet_Gas_Temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Outlet_Gas_Temp as text
%        str2double(get(hObject,'String')) returns contents of Outlet_Gas_Temp as a double


% --- Executes during object creation, after setting all properties.
function Outlet_Gas_Temp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Outlet_Gas_Temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Outlet_Liquid_Temp_Callback(hObject, eventdata, handles)
% hObject    handle to Outlet_Liquid_Temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Outlet_Liquid_Temp as text
%        str2double(get(hObject,'String')) returns contents of Outlet_Liquid_Temp as a double


% --- Executes during object creation, after setting all properties.
function Outlet_Liquid_Temp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Outlet_Liquid_Temp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
