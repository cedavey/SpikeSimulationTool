function varargout = SST(varargin)
% SST MATLAB code for SST.fig
%      SST, by itself, creates a new SST or raises the existing
%      singleton*.
%
%      H = SST returns the handle to a new SST or the handle to
%      the existing singleton*.
%
%      SST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SST.M with the given input arguments.
%
%      SST('Property','Value',...) creates a new SST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SST_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SST_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SST

% Last Modified by GUIDE v2.5 12-Nov-2019 13:52:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SST_OpeningFcn, ...
                   'gui_OutputFcn',  @SST_OutputFcn, ...
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


% --- Executes just before SST is made visible.
function SST_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SST (see VARARGIN)

% Choose default command line output for SST
handles.output = hObject;
handles = setSstTooltips(handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SST wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SST_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in runButton.
function runButton_Callback(hObject, eventdata, handles)
data = struct;
data.Naxons      = str2num(handles.nAxonsTextBox.String);
data.SNR         = str2num(handles.snrTextBox.String);
data.total_time  = str2num(handles.totalTimeTextBox.String);
data.fs          = str2num(handles.samplingRateTextBox.String);
data.sr          = randi(str2num(handles.maxSpikeRateTextBox.String),1,data.Naxons);
data.overlap     = handles.overlapCheckBox.Value;
data.rpt_temp    = handles.repeatTemplateCheckBox.Value;
data.has_drift   = handles.hasDriftCheckBox.Value;
data.has_noise   = handles.hasNoiseCheckBox.Value;
data.pre_noise   = handles.preNoiseCheckBox.Value;
data.do_filter   = handles.doFilterCheckBox.Value;
data.passband    = [str2num(handles.lowBandTextBox.String) str2num(handles.highBandTextBox.String)];
data.PLOT        = handles.plotCheckBox.Value;

SpikeSimulationTool(data);

% --- Executes during object creation, after setting all properties.
function nAxonsTextBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function samplingRateTextBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function totalTimeTextBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function maxSpikeRateTextBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function snrTextBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function lowBandTextBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function highBandTextBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
