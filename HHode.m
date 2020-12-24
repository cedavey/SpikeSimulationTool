%HH.m

function xdot= HHode(t, x, flag, Iapp)
% ODE variables
V=x(1); m=x(2); h=x(3); n=x(4);

% Set constant parameter (Can be changed for different cell types) 
VNa = 50;      gNa = 120; %Sodium equi voltage & conductance
VK  = -77;     gK  = 36;  %Potassium equi voltage & conductance
VL  = -54.4;   gL  = 0.3; %Leak current equi voltage & conductance
C   = 1;  %Capacitance of membrane

%taum= 1/ (         alpha(V)              +      beta(V)     )
taum = 1/ (0.1*(V+40)/(1-exp(-(V+40)/10)) + 4*exp(-(V+65)/18));
%minf=            alpha(V)             *taum
minf = (0.1*(V+40)/(1-exp(-(V+40)/10)))*taum;

tauh = 1/(0.07*exp(-(V+65)/20)+(1/(1+exp(-(V+35)/10))));
hinf = 0.07*exp(-(V+65)/20)*tauh;

taun = 1/(0.01*(V+55)/(1-exp(-(V+55)/10))+0.125*exp(-(V+65)/80));
ninf = (0.01*(V+55)/(1-exp(-(V+55)/10)))*taun;

%

xdot(1,1) = (Iapp(t) - gNa*m^3*h*(V-VNa) - gK*(V-VK)*n^4 - gL*(V-VL))/C; %dV/dt
xdot(2,1) = -(m - minf)/taum; %dm/dt
xdot(3,1) = -(h - hinf)/tauh; %dh/dt
xdot(4,1) = -(n - ninf)/taun; %dn/dt
% End HH.m