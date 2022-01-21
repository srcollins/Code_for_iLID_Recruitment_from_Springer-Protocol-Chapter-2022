%% Script exploring iLID_PDE_global_model

%%
kDdark=4.7;%4.7 for WT, 47 for Micro;
kDlight=0.13;%0.13 for WT, 0.8 for Micro;
t_range=0:0.1:10;

kRevert = 0.02;  %iLID inactivation rate
%kOffLit = 0.5;  %iLID-SspB disassociation rate in the dark state
kOffLit = 0.25;  %iLID-SspB disassociation rate in the dark state
kBind = kOffLit/kDlight;  %iLID-SspB association rate
kOffDark = kBind*kDdark;  %iLID-SspB association rate in the active/lit state
SspBTot = 0.1;  %total concentration of SspB
iLIDTot = 0.55;   %total concentration of iLID
D=25;   %sspB cytosolic diffusion rate
cellRadius = 5;

x_range=0:0.1:cellRadius;
shellFactors=x_range.^2 * 4 * pi;

p0=[kRevert kBind D SspBTot iLIDTot cellRadius]
%% Default parameters - ODE
pODE(1) = kRevert;  % iLID inactivation rate
pODE(3) = kOffLit;  % disassociation iLID(active)-SspB  (Kd dark estimate of 0.13 uM from https://www.pnas.org/content/115/10/E2238.long
pODE(2) = pODE(3)/kDlight;     % association iLID(active)-SspB  (t1/2 association estimate of 2.0 sec from https://www.pnas.org/content/115/10/E2238.long)
pODE(4) = pODE(2);     % association iLID(inactive)-SspB
pODE(5) = pODE(4)*kDdark;     % disassociation iLID(inactive)-SspB (Kd dark estimate of 4.7 uM from https://www.pnas.org/content/115/10/E2238.long
pODE(6) = SspBTot;    % total SspB concentration (microMolar)
pODE(7) = iLIDTot;    % total iLID concentration
inputStruct=makeInputStepFunStruct(0.1,0.1,1000);
[f,yAll] = iLID_ODE_model(pODE,t_range, inputStruct);
%%
[t,x,u]=iLID_PDE_global_model(t_range,p0,x_range);
% u: 1st dimension is time
% u: 2nd dimension is space
% u: 3rd dimension is component
%%
figure;
leg1=arrayfun(@(x) {[sprintf('%.1f',(x-1)/10) ' microns']},1:25:length(x));
plot(t,u(:,[1:25:length(x)],5)); hold on;
plot(t-0.1,yAll(:,5),'LineWidth',2,'Color','k');
legend([leg1 {'ODE'}]);
xlim([0 max(t)]);
%%
rad=1:10;
for i=1:length(rad)
    x_range=0:0.1:rad(i);
    shellFactors=x_range.^2 * 4 * pi;
    maxFactor(i)=max(shellFactors);
    maxFactorRel(i)=max(shellFactors)/sum(shellFactors);
    approxFactor(i)=(4*pi*rad(i)^2) / (4/3*pi*rad(i)^3);
end
plot(rad,1./maxFactorRel); hold on;
plot(rad,10*rad/3);

%% Titrate iLID and compare recruitment and relaxation kinetic profiles
%SspBTot=0.5;
%iLIDconcArr=[0.1 0.25 0.5 1];
SspBTot=1.5;
%iLIDconcArr=[0.1 0.25 0.5 1];
iLIDconcArr=[0.05 1.5 2.0 2.5];
cellRadius = 15;
x_range=0:0.1:cellRadius;
t_range=0:0.1:60;
f=nan(length(iLIDconcArr),length(t_range));

for i=1:length(iLIDconcArr)
    pin=[kRevert kBind D SspBTot iLIDconcArr(i) cellRadius];
    [t,x,u]=iLID_PDE_global_model(t_range,pin,x_range);
    f(i,:)=u(:,1,5);
end
%%
figure; 
colors=parula(length(iLIDconcArr));
subplot(1,2,1); hold on;
for i=1:length(iLIDconcArr)
    plot(t,f(i,:),'Color',colors(i,:),'LineWidth',2);
end
subplot(1,2,2);
hold on;
for i=1:length(iLIDconcArr)
    v=f(i,:);
    v=(v-min(v))/(max(v)-min(v));
    plot(t,v,'Color',colors(i,:),'LineWidth',2);
end

%%

figure;
colors=parula(length(iLIDconcArr));
hold on;
for i=1:length(iLIDconcArr)
    v=f(i,:);
    v=(v-min(v))/(max(v)-min(v));
    plot(t,v,'Color',colors(i,:),'LineWidth',3);
end

ylim([-0.1 1.1]);
xlabel('Time (sec)');
ylabel('Normalized cyto. intensity');

%%

