
function xdot= HHode(t, x, flag, Iapp, parameters)
% ODE variables
V=x(1); m=x(2); h=x(3); n=x(4);

% Calculates all the alphas and beta
nAlpha = 0.01*(V - (parameters.vRest + 10))/(1 - exp(-(V - (parameters.vRest + 10))/10));
nBeta = 0.125*exp(-(V - parameters.vRest)/80);
mAlpha = 0.1*(V - (parameters.vRest + 25))/(1 - exp(-(V - (parameters.vRest + 25))/10));
mBeta = 4*exp(-(V - parameters.vRest)/18);
hAlpha = 0.07*exp(-(V - parameters.vRest)/20);
hBeta = 1/(1 + exp(-(V - (parameters.vRest + 30))/10));

% Calculates tau and inf
tauN = 1/(nAlpha + nBeta);
infN = nAlpha * tauN;

tauM = 1/(mAlpha + mBeta);
infM = mAlpha * tauM;

tauH = 1/(hAlpha + hBeta);
infH = hAlpha * tauH;

% calculates the ODE
xdot(1,1) = (Iapp(t) - parameters.gNa * m^3 * h * (V - parameters.eNa) - parameters.gK * (V - parameters.eK) * n^4 - parameters.gLeak * (V - parameters.eLeak))/parameters.C;
xdot(2,1) = -(m - infM)/tauM;
xdot(3,1) = -(h - infH)/tauH;
xdot(4,1) = -(n - infN)/tauN;
end