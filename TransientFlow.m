function varargout = TransientFlow(varargin)
% TRANSIENTFLOW M-file for TransientFlow.fig
%      TRANSIENTFLOW, by itself, creates a new TRANSIENTFLOW or raises the existing
%      singleton*.
%
%      H = TRANSIENTFLOW returns the handle to a new TRANSIENTFLOW or the handle to
%      the existing singleton*.
%
%      TRANSIENTFLOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRANSIENTFLOW.M with the given input arguments.
%
%      TRANSIENTFLOW('Property','Value',...) creates a new TRANSIENTFLOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TransientFlow_OpeningFcn gets called.  An
%      unrecognized property name or i  nvalid value makes property application
%      stop.  All inputs are passed to TransientFlow_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TransientFlow

% Last Modified by GUIDE v2.5 11-Nov-2008 09:01:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @TransientFlow_OpeningFcn, ...
    'gui_OutputFcn',  @TransientFlow_OutputFcn, ...
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


% --- Executes just before TransientFlow is made visible.
function TransientFlow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TransientFlow (see VARARGIN)

% Choose default command line output for TransientFlow
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TransientFlow wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.uitable1,'Data',infiltration(1));

data=get(handles.uitable1,'Data');
data(isnan(data(:,1)),:)=[];
data(isnan(data(:,2)),:)=[];
axes(handles.axes1)
plot(data(:,1),data(:,2),'-b','linewidth',1.5)
xlabel('Time (d)');
ylabel('Infiltration (mm/d)');
box on; grid on;

uitable2_CellEditCallback(hObject, eventdata, handles)

set(handles.text3,'String','Hydrostatic');
set(handles.text5,'String','');

% --- Outputs from this function are returned to the command line.
function varargout = TransientFlow_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

data=get(handles.uitable1,'Data');
data(isnan(data(:,1)),:)=[];
data(isnan(data(:,2)),:)=[];
axes(handles.axes1)
hold off
plot(data(:,1),data(:,2),'-b','linewidth',1.5)
xlabel('Time (d)');
ylabel('Infiltration (mm/d)');
box on; grid on;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc

data=get(handles.uitable1,'Data');
data(isnan(data(:,1)),:)=[];
data(isnan(data(:,2)),:)=[];
tB=data(:,1);
qT=-data(:,2)/1000;

data=get(handles.uitable2,'Data');
SoilPars.thetaR=data(1);
SoilPars.thetaS=data(2);
% SoilPars.alpha=data(3);
% SoilPars.n=data(4);
SoilPars.psiS=data(3);
SoilPars.lambda=data(4);

SoilPars.Ks=data(5);
SoilPars.neta=data(6);
zN=data(7);
t=0:data(9):data(8);

zWT=0;

% if PlotHydProps
%     CheckHydProps(SoilPars);
% end

% Finite difference grid space step:
dz=zN/100;

psiT=[];
qB=[];        % Lower boundary flux [m/d], or
psiB=zWT;     %                               lower boundary psi [m]

% Initial conditions:
zIC=[0 zN];
psiIC=[zWT zWT-zN];
% psiIC=[0 0];

% Finite difference grid:
z=dz/2:dz:zN-dz/2;
n=numel(z);

% Initial condition:
psi0=interp1(zIC,psiIC,z);

% Solve Richards equation using an ODE solver:
JPat=spdiags(ones(n,3),[-1 0 1],n,n);
options=odeset('JPattern',JPat,'MaxStep',1);
wait=waitbar(0,'Flow model is solving...');
[t,psi]=ode15s(@Richards,t,psi0,options,[t(1) t(end)],dz,n,SoilPars,psiB,psiT,qB,qT,tB);
close(wait)

% Get water content:
theta=thetaFun(psi,SoilPars);

% Get hydraulic head:
h=psi+repmat(z,numel(t),1);

hB=h(:,1);
Kb=KFun(hB/2,SoilPars);
rchg=1000*(-Kb.*(0-hB)/(dz/2));

dt=t(2)-t(1);
% S=1000*[0; diff(sum(theta*dz,2))/dt];
S=1000*sum(theta*dz,2);
S=S-S(1);
% S=[0

clear data;
data.psi=psi;
data.h=h;
data.theta=theta;
data.t=t;
data.z=z;
data.SoilPars=SoilPars;
data.tB=tB;
data.qT=qT;
data.rchg=rchg;
data.S=S;
save data data

i=1;
set(handles.edit1,'String',num2str(i));
UpdatePlot(i,hObject, eventdata, handles)

% set(handles.text5,'String','Recharge');
% for i=1:numel(t);
%     set(handles.text3,'String',['t =' num2str(t(i)) ' d']);
% 
%     axes(handles.axes2)
%     hold off
%     plot(psi(i,:),z,'-b','linewidth',1.5)
%     hold on
%     plot(z,z,'-g','linewidth',1.5)
%     plot(h(i,:),z,'-r','linewidth',1.5)
%     legend('Pressure head','Elevation head','Hydraulic head',4)
%     xlabel('Head (m)');
%     ylabel('Elevation (m)');
%     box on; grid on;
%     xmin=floor(max(min([min(psi(:)) min(h(:)) min(z)]),-10));
%     xmax=ceil(min(max([max(psi(:)) max(h(:)) max(z)]),10));
%     xlim([xmin xmax]);
%     
%    
%     axes(handles.axes3)
%     hold off
%     plot(theta(i,:),z,'-m','linewidth',1.5)
%     hold on
%     plot(SoilPars.thetaR*[1 1],ylim,':r')
%     plot(SoilPars.thetaS*[1 1],ylim,':r')
%     xlabel('Water content (-)');
%     ylabel('Elevation (m)');
%     box on; grid on;
%     xlim([0 ceil(SoilPars.thetaS*10+1)/10]);
%     
%     axes(handles.axes1)
%     v=ylim;
%     hold off
%     plot(tB,-qT*1000,'-b',t,rchg,'-r','linewidth',1.5)
%     hold on
%     area([0 t(i)],[100 100],-100,'facecolor',[0.8 0.8 0.8],'linestyle','none');
%     plot(tB,-qT*1000,'-b',t,rchg,'-r','linewidth',1.5)
%     xlabel('Time (d)');
%     ylabel('Flux (mm/d)');
%     box on; grid on;
% 	ylim(v)
% 	xlim([0 t(end)]);
% %     legend('Infiltration','Recharge',1);
%     
%     pause(0.1);
% end


function [dpsidt]=Richards(t,psi,tspan,dz,n,SoilPars,psiB,psiT,qB,qT,tB)
%Solves Richards equation in one dimension.
%(Richards, 1931, Physics, 1, 318-333)
%
%C [L-1] - Specific capacity
%dz [L] - Space step
%psi [L] - Pressure head
%psiB [L] - Pressure head at base
%psiT [L] - Pressure head at top
%K [LT-1] - Hydraulic conductivity
%KIn [LT-1] - Hydraulic conductivity at base of block
%KOut [LT-1] - Hydraulic conductivity at top of block
%Se [-] - Effective saturation
%z [L] - Depth
%qB [LT-1] - Flux at base
%qT [LT-1] - Flux at top
%q [LT-1] - Flux
%qIn [LT-1] - Flux into base of block
%qOut [LT-1] - Flux out of top of block

%Define the i index
i=1:n-1;
%Calculate hydraulic conductivity
Kn=KFun(psi,SoilPars);
K(i+1,1)=(Kn(i+1)+Kn(i,1))/2;
%Calculate hydraulic head gradients:
dpsidz(i+1,1)=(psi(i+1)-psi(i))/dz+1;
%Calculate flows:
q(i+1,1)=-K(i+1).*dpsidz(i+1);

%Lower boundary:
if ~isnan(qB)
    q(1)=qB;
else
    KB=(Kn(1)+KFun(psiB,SoilPars))/2;
    dpsidz(1)=(psi(1)-psiB)/dz*2+1;
    q(1)=-KB*dpsidz(1);
end

%Upper boundary:
if ~isnan(qT)
    if numel(qT)>1;
        qT=interp1(tB,qT,t);
    end
    q(n+1)=qT;
else
    KT=(Kn(n)+KFun(psiT,SoilPars))/2;
    dpsidz(n+1)=(psiT-psi(n))/dz*2+1;
    q(n+1)=-KT*dpsidz(n+1);
end

%Apply continuity
dthetadt=-diff(q)/dz;

%Calculate specific capacity
C=CFun(psi,SoilPars);

%Calculate pressure gradient:
dpsidt=dthetadt./C;

waitbar((t-tspan(1))/(tspan(2)-tspan(1)));

function C=CFun(psi,SoilPars)
%Calculates specific capacity, C as a function of head, psi
%(van GenucpsiTen, 1980, Soil Sci. Soc. Am. J. 44, 892–898)
Ss=1e-3;
thetaR=SoilPars.thetaR;
thetaS=SoilPars.thetaS;
% alpha=SoilPars.alpha;
% n=SoilPars.n;
% m=1-1/n;;
% Se=(1+abs(alpha.*psi).^n).^-m;
% Se(psi>=0)=1;
% dSedh=alpha.*m./(1-m).*Se.^(1./m).*(1-Se.^(1./m)).^m;
% dSedh(psi>=0)=0;
psiS=SoilPars.psiS;
lambda=SoilPars.lambda;
Se=(psiS./psi).^lambda;
Se(psi>=psiS)=1;
dSedh=-lambda./psi.*Se;
dSedh(psi>=psiS)=0;

C=Se.*Ss+(thetaS-thetaR).*dSedh;

%**************************************************************************

function K=KFun(psi,SoilPars)
%Calculates hydraulic conductivity, K as a
%function of head, effective saturation, Se

%Mualem relationship:
thetaR=SoilPars.thetaR;
thetaS=SoilPars.thetaS;
Ks=SoilPars.Ks;
neta=SoilPars.neta;
% alpha=SoilPars.alpha;
% n=SoilPars.n;
% m=1-1/n;
% Se=(1+abs(alpha.*psi).^n).^-m;
% kr=Se.^neta.*(1-(1-Se.^(1./m)).^m).^2;

psiS=SoilPars.psiS;
lambda=SoilPars.lambda;
Se=(psiS./psi).^lambda;
Se(psi>=psiS)=1;
kr=Se.^(neta+2+2/lambda);

K=Ks.*kr;

%**************************************************************************

function [theta]=thetaFun(psi,SoilPars)
%Calculates moisture content, theta and effective saturation, Se as a
%function of head, psi
%(van GenucpsiTen, 1980, Soil Sci. Soc. Am. J. 44, 892–898)
thetaR=SoilPars.thetaR;
thetaS=SoilPars.thetaS;

% alpha=SoilPars.alpha;
% n=SoilPars.n;
% m=1-1/n;
% Se=(1+abs(alpha.*psi).^n).^-m;

psiS=SoilPars.psiS;
lambda=SoilPars.lambda;
Se=(psiS./psi).^lambda;
Se(psi>=psiS)=1;

theta=(thetaS-thetaR).*Se+thetaR;

%**************************************************************************

function CheckHydProps(SoilPars);

psi=-logspace(-2,2);
th=thetaFun(psi,SoilPars);
C=CFun(psi,SoilPars);
K=KFun(psi,SoilPars)./SoilPars.Ks;

figure(2)
subplot(3,1,1)
semilogx(psi,th,'linewidth',1.5);
ylabel('\theta (-)')
title('Hydraulic properties','fontsize',14)
subplot(3,1,2);
semilogx(psi,C,'linewidth',1.5);
ylabel('C (m^{-1})')
subplot(3,1,3);
semilogx(psi,K,'linewidth',1.5);
ylabel('K (m/d)')
xlabel('\psi (m)')

set(gcf,'position',[850 100 400 600]);

% --- Executes when entered data in editable cell(s) in uitable2.
function uitable2_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

clc

data=get(handles.uitable2,'Data');
SoilPars.thetaR=data(1);
SoilPars.thetaS=data(2);
% SoilPars.alpha=data(3);
% SoilPars.n=data(4);
SoilPars.psiS=data(3);
SoilPars.lambda=data(4);

SoilPars.Ks=data(5);
SoilPars.neta=data(6);
zN=data(7);
t=0:data(9):data(8);

z=linspace(0,zN,100);
psi=-z;
h=zeros(size(z));

theta=thetaFun(psi,SoilPars);
i=1;
axes(handles.axes2)
hold off
plot(psi(i,:),z,'-b','linewidth',1.5)
hold on
plot(z,z,'-g','linewidth',1.5)
plot(h(i,:),z,'-r','linewidth',1.5)
xlabel('Head (m)');
ylabel('Elevation (m)');
box on; grid on;
xmin=floor(min([min(psi(:)) min(h(:)) min(z)]));
xmax=ceil(max([max(psi(:)) max(h(:)) max(z)]));
xlim([xmin xmax]);
legend('Pressure head','Elevation head','Hydraulic head')

axes(handles.axes3)
hold off
plot(theta(i,:),z,'-m','linewidth',1.5)
hold on
plot(SoilPars.thetaR*[1 1],ylim,':r')
plot(SoilPars.thetaS*[1 1],ylim,':r')
xlabel('Water content (-)');
ylabel('Elevation (m)');
box on; grid on;
xlim([0 ceil(SoilPars.thetaS*10+1)/10]);

set(handles.text3,'String','Hydrostatic');

function data=infiltration(r)

data=[         0         0         0         0         0         0         0         0         0         0
    0.9990         0    0.9990         0    0.9990         0    0.9990         0    1.0000   12.7211
    1.0000   10.0000    1.0000   40.0000    1.0000   -0.1000    1.0000   20.0000    2.0000    5.4297
   10.0000   10.0000    1.9990   40.0000   100.000   -0.1000    1.9990   20.0000    3.0000   -2.4580
       NaN       NaN    2.0000         0       NaN       NaN    2.0000   -2.0000    4.0000    1.2885
       NaN       NaN   10.0000         0       NaN       NaN    2.9990   -2.0000    5.0000   -0.3079
       NaN       NaN       NaN       NaN       NaN       NaN    3.0000         0    6.0000   16.7523
       NaN       NaN       NaN       NaN       NaN       NaN   20.0000         0    7.0000    1.9220
       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN    8.0000    8.7404
       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN    9.0000    4.0379
       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN   10.0000    9.4431
       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN   11.0000    1.1758
       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN   12.0000   -2.0113
       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN   13.0000   -2.1119
       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN   14.0000    7.4410
       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN   15.0000    4.1131
       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN   16.0000    4.0197
       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN   17.0000   12.0966
       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN   18.0000    6.4579
       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN   19.0000    5.9891
       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN   20.0000   12.9385];

data=data(:,2*r+[-1 0]);
if r==5
    data(:,2)=randn(21,1)*5+5;
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
r=get(handles.popupmenu1,'Value')
set(handles.uitable1,'Data',infiltration(r));
uitable1_CellEditCallback(hObject, eventdata, handles)

if r==1
    data=get(handles.uitable2,'Data');
    data(8)=10;
    data(9)=0.1;
    set(handles.uitable2,'Data',data);
elseif r==2
    data=get(handles.uitable2,'Data');
    data(8)=10;
    data(9)=0.1;
    set(handles.uitable2,'Data',data);
elseif r==3
    data=get(handles.uitable2,'Data');
    data(8)=100;
    data(9)=1;
    set(handles.uitable2,'Data',data);
elseif r==4
    data=get(handles.uitable2,'Data');
    data(8)=20;
    data(9)=0.2;
    set(handles.uitable2,'Data',data);
elseif r==5
    data=get(handles.uitable2,'Data');
    data(8)=20;
    data(9)=0.25;
    set(handles.uitable2,'Data',data);
end
% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', {'S.S. Infiltration', 'Rainfall pulse', 'Continuous Evaporation', 'Rain - Evap', 'Random'});



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton5,'Enable','Off');
i=max(1,str2num(get(handles.edit1,'String'))-1);
set(handles.edit1,'String',num2str(i));
UpdatePlot(i,hObject, eventdata, handles)
set(handles.pushbutton5,'Enable','On');
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton6,'Enable','Off');
i=str2num(get(handles.edit1,'String'))+1;
set(handles.edit1,'String',num2str(i));
UpdatePlot(i,hObject, eventdata, handles)
set(handles.pushbutton6,'Enable','On');
% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton7,'Enable','Off');
i=max(1,str2num(get(handles.edit1,'String'))-10);
set(handles.edit1,'String',num2str(i));
UpdatePlot(i,hObject, eventdata, handles)
set(handles.pushbutton7,'Enable','On');
% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton8,'Enable','Off');
i=str2num(get(handles.edit1,'String'))+10;
set(handles.edit1,'String',num2str(i));
UpdatePlot(i,hObject, eventdata, handles)
set(handles.pushbutton8,'Enable','On');

function UpdatePlot(i,hObject, eventdata, handles)

load data;
t=data.t;
psi=data.psi;
theta=data.theta;
S=data.S;
h=data.h;
z=data.z;
SoilPars=data.SoilPars;
tB=data.tB;
qT=data.qT;
rchg=data.rchg;

if i>numel(t)
    i=numel(t);
end

set(handles.text5,'String','Drainage');
set(handles.text3,'String',['t =' num2str(t(i)) ' d']);

axes(handles.axes2)
hold off
plot(psi(i,:),z,'-b','linewidth',1.5)
hold on
plot(z,z,'-g','linewidth',1.5)
plot(h(i,:),z,'-r','linewidth',1.5)
legend('Pressure head','Elevation head','Hydraulic head')
xlabel('Head (m)');
ylabel('Elevation (m)');
box on; grid on;
xmin=floor(max(min([min(psi(:)) min(h(:)) min(z)]),-10));
xmax=ceil(min(max([max(psi(:)) max(h(:)) max(z)]),10));
xlim([xmin xmax]);


axes(handles.axes3)
hold off
plot(theta(i,:),z,'-m','linewidth',1.5)
hold on
plot(SoilPars.thetaR*[1 1],ylim,':r')
plot(SoilPars.thetaS*[1 1],ylim,':r')
xlabel('Water content (-)');
ylabel('Elevation (m)');
box on; grid on;
xlim([0 ceil(SoilPars.thetaS*10+1)/10]);

axes(handles.axes1)
% v=ylim;
hold off
plot(tB,-qT*1000,'-b',t,rchg,'-r','linewidth',1.5)
hold on
area([0 t(i)],[100 100],-100,'facecolor',[0.8 0.8 0.8],'linestyle','none');
plot(tB,-qT*1000,'-b',t,rchg,'-r',t,S,'-g','linewidth',1.5)
xlabel('Time (d)');
ylabel('Flux (mm/d)');
box on; grid on;
ylim([min([S;(-qT*1000);rchg]) max([S;(-qT*1000);rchg])])
% ylim(v)
xlim([0 t(end)]);





