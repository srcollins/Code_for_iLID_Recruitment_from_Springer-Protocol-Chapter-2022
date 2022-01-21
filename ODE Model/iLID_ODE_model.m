function [f,yAll] = iLID_ODE_model(p,t,inputStruct,relTol)
% Basic iLID model

% Species
% 1 = Free iLID inactive
% 2 = Free iLID active 
% 3 = iLID-SspB (iLID inactive)
% 4 = iLID-SspB (iLID active)
% 5 = Free SspB

% Rate constants
if isempty(p)
    kDdark=4.7;%4.7 for WT, 47 for Micro;
    kDlight=0.13;%0.13 for WT, 0.8 for Micro;
    p(1) = 0.02;  % iLID inactivation rate
    p(2) = 0.5/kDlight;     % association iLID(active)-SspB  (t1/2 association estimate of 2.0 sec from https://www.pnas.org/content/115/10/E2238.long)
    p(3) = 0.5;  % disassociation iLID(active)-SspB  (Kd dark estimate of 0.13 uM from https://www.pnas.org/content/115/10/E2238.long
    p(4) = p(2);     % association iLID(inactive)-SspB
    p(5) = p(4)*kDdark;     % disassociation iLID(inactive)-SspB (Kd dark estimate of 4.7 uM from https://www.pnas.org/content/115/10/E2238.long
    p(6) = 0.5;    % total SspB concentration (microMolar)
    p(7) = 0.1;    % total iLID concentration
end
% Initial conditions
iLIDtot=p(7);
b=iLIDtot+p(6)+p(5)/p(4);
y0(2)=0;
y0(3)=(b - sqrt(b^2 - 4*iLIDtot*p(6)))/2;
y0(1)=iLIDtot-y0(3);
y0(4)=0;
y0(5)=p(6)-y0(3);

for i=1:length(inputStruct.startTime)
    dy{i} = @(t,y) [p(1)*y(2) + p(5)*y(3) - inputStruct.inputStrength(i)*y(1) - p(4)*y(1)*y(5);     % 1 = Free iLID inactive
             inputStruct.inputStrength(i)*y(1) + p(3)*y(4) - p(1)*y(2) - p(2)*y(2)*y(5);            % 2 = Free iLID active 
             p(4)*y(1)*y(5) + p(1)*y(4) - inputStruct.inputStrength(i)*y(3) - p(5)*y(3);            % 3 = iLID-SspB (iLID inactive)
             inputStruct.inputStrength(i)*y(3) + p(2)*y(2)*y(5) - p(1)*y(4) - p(3)*y(4);            % 4 = iLID-SspB (iLID active)
             p(3)*y(4) + p(5)*y(3) - p(2)*y(2)*y(5) - p(4)*y(1)*y(5)];                              % 5 = Free SspB
end

yAll=nan(length(t),length(y0));
tAll=union(t,inputStruct.startTime);
if nargin<4
    t0=0;
    for i=1:length(inputStruct.startTime)
        tF=min([inputStruct.startTime((i+1):end) max(tAll)]);
        tSegment=tAll(tAll>=t0 & tAll<=tF);
        [~,ySegment] = ode45(dy{i},tSegment,y0);
        y0=ySegment(end,:);
        t0=tF;
        [~,i1,i2]=intersect(t,tSegment);
        yAll(i1,:)=ySegment(i2,:);
    end
else % Relative Tolerance for ode45 was provided as an extra input
    options=odeset('RelTol',relTol,'AbsTol',1e-8);
    t0=0;
    for i=1:length(inputStruct.startTime)
        tF=min([inputStruct.startTime((i+1):end) max(tAll)]);
        tSegment=tAll(tAll>=t0 & tAll<=tF);
        [~,ySegment] = ode45(dy{i},tSegment,y0,options);      % Call ode45 multiple times for receptor on and off
        y0=ySegment(end,:);
        t0=tF;
        [~,i1,i2]=intersect(t,tSegment);
        yAll(i1,:)=ySegment(i2,:);
    end
end
f=yAll(:,3)+yAll(:,4);
