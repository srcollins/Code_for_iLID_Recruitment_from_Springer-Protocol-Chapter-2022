function [t,x,u]=iLID_PDE_global_model(t_range,p,x_range,LightInput)
%this program simulates a model of iLID binding to sspB in a spherical
%cell, considering local diffusion of sspB as a potential limiting factor

% Species
% 1 = Free iLID inactive
% 2 = Free iLID active 
% 3 = iLID-SspB (iLID inactive)
% 4 = iLID-SspB (iLID active)
% 5 = Free SspB

%Einstein equation estimates that the density should propagate a distance
%of sqrt(4*D*t) over time t

global D;   %cytoplasmic diffusion coefficient of sspB
global kRevert;  %iLID inactivation rate
global kBind;  %iLID-SspB association rate
global kOffDark;  %iLID-SspB association rate in the dark state
global kOffLit;  %iLID-SspB association rate in the active/lit state
global SspBTot;  %total concentration of SspB
global iLIDTot;   %total concentration of iLID
global kDdark;  %iLID-SspB Kd in the dark state
global cellRadius; 

%Default parameter values
D=[0 0 0 0 25]';         % Estimates from the literature: 1 for PIP2, 3 for PIP3, 0.1 for G-proteins
                        % .001 for a slow TM receptor, .05 for a mobile TM
                        % receptor, .0045 to .003 for LDL Receptor, .05 for
                        % Fc epsilon receptor, 25 for GFP in a cell

kDdark=4.7;%4.7 for WT, 47 for Micro;
kDlight=0.13;%0.13 for WT, 0.8 for Micro;

D = 25;   %vector of diffusion coefficients
kRevert = 0.02;  %iLID inactivation rate
kBind = 0.5/kDlight;  %iLID-SspB association rate
kOffLit = 0.5;  %iLID-SspB association rate in the dark state
kOffDark = kBind*kDdark;  %iLID-SspB association rate in the active/lit state
SspBTot = 0.5;  %total concentration of SspB
iLIDTot = 0.1;   %total concentration of iLID
cellRadius = 10;

if nargin>1 & ~isempty(p)
    % Parameters:
    % p(1) = iLID reversion rate
    % p(2) = iLID-SspB association rate
    % p(3) = sspB cytosolic diffusion rate
    % p(4) = SspBTot
    % p(5) = iLIDTot
    D = [zeros(1,4) p(3)]';   %vector of diffusion coefficients
    kRevert = p(1);  %iLID inactivation rate
    kBind = p(2);  %iLID-SspB association rate
    kOffDark = kBind*kDdark;  %iLID-SspB association rate in the dark state
    kOffLit = kBind*kDlight;  %iLID-SspB association rate in the active/lit state
    SspBTot = p(4);  %total concentration of SspB
    iLIDTot = p(5);   %total concentration of iLID
    cellRadius = p(6);  % cell radius in microns
end
if nargin<3
    x_range=0:0.1:cellRadius; %distance from cell centroid to periphery in microns
else
    if max(x_range) ~= cellRadius
        fprintf(' **** xRange does not match cell radius ***\n\n');
    end
end

% Set Receptor input spatial pattern
if nargin<4
    R=@ReceptorInputFunction;
else
    R=ReceptorInput;
end


m = 2;         %symmetry, 0=slab, 1=cylindrical, 2=spherical
x = x_range;   %radius of the cell, in microns
t = t_range;

u = pdepe(m,@diff1Dpde,@diff1Dic,@diff1Dbc,x,t);


% --------------------------------------------------------------------------

function [c,f,s] = diff1Dpde(x,t,u,DuDx)
global D;   %vector of diffusion coefficients
global kRevert;  %iLID inactivation rate
global kBind;  %iLID-SspB association rate
global kOffDark;  %iLID-SspB association rate in the dark state
global kOffLit;  %iLID-SspB association rate in the active/lit state
global SspBTot;  %total concentration of SspB
global iLIDTot;   %total concentration of iLID
global kDdark;  %iLID-SspB Kd in the dark state

c = ones(5,1);
f = D .* DuDx;  % D should be a column vector

s(1,1) = kRevert*u(2) + kOffDark*u(3) - kBind*u(1)*u(5);
s(2,1) = kOffLit*u(4) - kRevert*u(2) - kBind*u(2)*u(5);
s(3,1) = kBind*u(1)*u(5) + kRevert*u(4) - kOffDark*u(3);
s(4,1) = kBind*u(2)*u(5) - kRevert*u(4) - kOffLit*u(4);
s(5,1) = kOffDark*u(3) + kOffLit*u(4) - kBind*u(1)*u(5) - kBind*u(2)*u(5);

% --------------------------------------------------------------------------

function u0 = diff1Dic(x)
% Defined to start with all iLID just activated
% Determine initial conditions
global iLIDTot;
global SspBTot;
global kDdark;
global cellRadius;

b=iLIDTot+SspBTot+kDdark;
basalBinding=(b - sqrt(b^2 - 4*iLIDTot*SspBTot))/2;
scaleFactor=2*10*cellRadius/3; %adjusts iLID conc at the membrane based on cell volume

u0(1,1) = 0; 
u0(2,1) = (iLIDTot - basalBinding) * (x==cellRadius) * scaleFactor;
u0(3,1) = 0;
u0(4,1) = basalBinding * (x==cellRadius) * scaleFactor;
u0(5,1) = SspBTot - basalBinding;
    
% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = diff1Dbc(xl,ul,xr,ur,t)
pl = zeros(5,1);
ql = ones(5,1);
pr = zeros(5,1);
qr = ones(5,1);

% --------------------------------------------------------------------------
function val=ReceptorInputFunction(x)
val=0;
if x>1 & x<=2
    val=1;
end
