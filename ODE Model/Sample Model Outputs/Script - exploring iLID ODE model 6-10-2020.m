%% Script for exploring the iLID ODE model

%% Default parameters
kDdark=4.7;%4.7 for WT, 47 for Micro;
kDlight=0.13;%0.13 for WT, 0.8 for Micro;
p0(1) = 0.02;  % iLID inactivation rate
p0(3) = 0.5;  % disassociation iLID(active)-SspB  (Kd dark estimate of 0.13 uM from https://www.pnas.org/content/115/10/E2238.long
p0(2) = p0(3)/kDlight;     % association iLID(active)-SspB  (t1/2 association estimate of 2.0 sec from https://www.pnas.org/content/115/10/E2238.long)
p0(4) = p0(2);     % association iLID(inactive)-SspB
p0(5) = p0(4)*kDdark;     % disassociation iLID(inactive)-SspB (Kd dark estimate of 4.7 uM from https://www.pnas.org/content/115/10/E2238.long
p0(6) = 0.25;    % total SspB concentration (microMolar)
p0(7) = 0.1;    % total iLID concentration
%%
t=0:0.2:100;
inputStruct=makeInputStepFunStruct(1,0.1,1000);
[f,yAll] = iLID_ODE_model(p0,t, inputStruct);
plot(t,yAll);
%% Titrate iLID conc
iLIDconc=10.^[-2:0.2:0.4];
inputStruct=makeInputStepFunStruct(1,0.1,1000);
t=0:0.2:100;
figure(1); hold on; figure(2); hold on;
fracRec=nan(length(iLIDconc),1);
tHalf=nan(length(iLIDconc),1);
for i=1:length(iLIDconc)
    pin=p0; pin(7)=iLIDconc(i);
    [f,yAll] = iLID_ODE_model(pin,t, inputStruct);
    fracRec(i)=1-min(yAll(:,5)/yAll(1,5));  % Minimum cytoplasmic concentration
    tHalf(i)=t(find(yAll(:,5)/yAll(1,5) < fracRec(i)+0.5*(1-fracRec(i)),1,'last'));
    figure(1); plot(t,yAll(:,5)/yAll(1,5))
    figure(2); plot(t,(yAll(:,5)-min(yAll(:,5)))/(yAll(1,5)-min(yAll(:,5))));
end
figure(3); scatter(fracRec,tHalf);

%% Titrate both concentrations
iLIDconc=10.^[-2:0.05:1];
sspbConc=10.^[-2:0.05:1];
t=0:0.2:200;
inputStruct=makeInputStepFunStruct(1,0.1,1000);
% structure for holding the simulation data
sim.t=t;
sim.pin=cell(length(iLIDconc),length(sspbConc));
sim.inputStruct=inputStruct;
sim.iLIDconc=iLIDconc;
sim.sspbConc=sspbConc;
sim.output=cell(length(iLIDconc),length(sspbConc));
sim.fracRec=nan(length(iLIDconc),length(sspbConc));
sim.absRec=nan(length(iLIDconc),length(sspbConc));
sim.basalRec=nan(length(iLIDconc),length(sspbConc));
sim.basalFracRec=nan(length(iLIDconc),length(sspbConc));
sim.tHalf=nan(length(iLIDconc),length(sspbConc));

for i=1:length(iLIDconc)
    for j=1:length(sspbConc)
        pin=p0; pin(7)=iLIDconc(i);
        pin(6)=sspbConc(j);
        sim.pin{i,j}=pin;
        [f,yAll] = iLID_ODE_model(pin,sim.t, inputStruct);
        sim.output{i,j}=yAll;
        sim.fracRec(i,j)=1-min(yAll(:,5)/sspbConc(j));  % Minimum cytoplasmic concentration
        sim.absRec(i,j)=max(yAll(:,3)+yAll(:,4));
        sim.basalRec(i,j)=min(yAll(:,3)+yAll(:,4));
        sim.basalFracRec(i,j)=sim.basalRec(i,j)/sspbConc(j);
        sim.tHalf(i,j)=t(find(yAll(:,5) < yAll(1,5) - 0.5*(sim.absRec(i,j)-sim.basalRec(i,j)),1,'last'));
    end
end
%% save results
save('G:\Shared drives\Collins lab\general-matlab\modeling\iLID\ODE global recruitment modeling data 06-19-2020.mat','sim');
%% load results
load('G:\Shared drives\Collins lab\general-matlab\modeling\iLID\ODE global recruitment modeling data 06-19-2020.mat','sim');
%%
tickInd=1:20:length(sspbConc);
figure;
subplot(2,3,1);
imagesc(sim.fracRec); colorbar;
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
xlabels=arrayfun(@(x) {num2str(x)},sspbConc(tickInd));
ylabels=arrayfun(@(x) {num2str(x)},iLIDconc(tickInd));
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Max Fraction Recruited')

subplot(2,3,2);
imagesc(log2(sim.tHalf)); colorbar;
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 tHalf for Dissociation');

subplot(2,3,3);
imagesc(log2(sim.absRec-sim.basalRec),log2([5e-4 5])); 
colorbar('Ticks',log2([1e-3 1e-2 0.1 1]),'TickLabels',[1e-3 1e-2 0.1 1]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Max Absolute Recruitment')

subplot(2,3,4);
imagesc(log2(sim.basalRec),log2([5e-4 5])); 
colorbar('Ticks',log2([1e-3 1e-2 0.1 1]),'TickLabels',[1e-3 1e-2 0.1 1]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Basal Recruitment');

subplot(2,3,5);
imagesc(sim.basalFracRec); colorbar;
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Basal Fraction of sspB Recruited');

subplot(2,3,6);
imagesc(log2(sim.absRec./sim.basalRec),log2([1 25]));
colorbar('Ticks',log2([1 2 5 10 20]),'TickLabels',[1 2 5 10 20]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Fold Enrichment by Recruitment');

%%
ilidConcMat=repmat(iLIDconc',[1 53]);
sspBconcMat=repmat(sspbConc,[53 1]);
figure;
scatter(log2(fracRec(:)./sspBconcMat(:)),tHalf(:));
ind1=sspBconcMat<0.25;
scatter(log2(fracRec(ind1)./sspBconcMat(ind1)),tHalf(ind1));
%scatter(log2(ilidConcMat(:)),fracRec(:));

%% plot time course graphs
for i=1:10:51
    plot(sim.t,sum(sim.output{i,i}(:,3:4),2));
    hold on
end
xlim([0 10]);
