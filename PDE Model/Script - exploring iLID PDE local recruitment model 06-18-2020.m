%% Script exploring iLID_PDE_model (local recruitment)

%% -----------------------------------------------------------------------
% Modeling / Fitting experimental data

%% Load experimental data
load('G:\Shared drives\Collins lab\general-matlab\modeling\iLID\experimental data\yOutFinal.mat','yOutFinal');
% 1 = CAAX
% 2 = Lyn
% 3 = Stargazin
% 4 = ADRB2
%%
x1=.2195*(0:100);  % x coordinate values for experimental data
profile1=nanmean(yOutFinal{3}{6});  % Use the stargazin data at t=1 to estimate the illumination pattern
profile1=(profile1 - 1)/(profile1(1)-1);
lightFun0=@(p,x) normpdf(x,0,p(1)) + p(2)*normpdf(x,0,p(3));
lightFun=@(p,x) lightFun0(p,x)/lightFun0(p,0);  % Normalize the function to have a max of 1
pLight=nlinfit(x1,profile1,@(p,x) lightFun(p,x), [3 0.2 10])

figure; hold on; plot(x1,profile1,'.'); plot(x1, lightFun(pLight,x1));
f=lightFun(pLight,x1);
ssr=sum((profile1-f).^2)

%% -----------------------------------------------------------------------
% Modeling the reaction diffusion system
kDdark=4.7;%4.7 for WT, 47 for Micro;
kDlight=0.13;%0.13 for WT, 0.8 for Micro;
t_range=0:0.1:10;

kRevert = 0.02;  %iLID inactivation rate
kOffLit = 0.5;  %iLID-SspB association rate in the dark state
%kOffLit = 0.75;  %iLID-SspB association rate in the dark state
kBind = kOffLit/kDlight;  %iLID-SspB association rate
kOffDark = kBind*kDdark;  %iLID-SspB association rate in the active/lit state
SspBTot = 0.5;  %total concentration of SspB
iLIDTot = 0.3;   %total concentration of iLID
D=0.1;   %iLID membrane lateral diffusion rate
cellRadius = 20;

x_range=0:0.1:cellRadius;

p0=[kRevert kBind D SspBTot iLIDTot cellRadius];

%%
%inputFun0=@(x) normpdf(x,0,3.3) + 0.2*normpdf(x,0,6);  % Earlier "eyeballed" estimate of the illumination pattern
%inputFun=@(x) inputFun0(x)/inputFun0(0);
inputFunFromFit=@(x) lightFun(pLight,x);
[t,x,u]=iLID_PDE_model(t_range,p0,x_range,inputFunFromFit);
% u: 1st dimension is time
% u: 2nd dimension is space
% u: 3rd dimension is component

%%
figure('Position',[200 400 1200 400]); hold on;
cmap=summer(ceil(length(t)/10));
cytoFrac=.075;
for i=1:10:length(t)
    if i==1
        yvals0=sum(u(i,:,3:4),3);
    end
    subplot(1,2,1);
    plot(x,sum(u(i,:,3:4),3),'Color',cmap(floor(i/10)+1,:),'LineWidth',2); hold on;
    if i==11 | i==101
        yvals=sum(u(i,:,3:4),3)-yvals0;
        xval=interp1(yvals,x,0.5*max(yvals))
        plot([xval xval],[0 0.5*max(yvals)+max(yvals0)],'k');
    end
    title('Recruitment (micromolar)')
    
    %Subpanel to approximate experimental results where fluorescence from
    %cytoplasm will contribute to the observed signal
    subplot(1,2,2);
    plot(x,(sum(u(i,:,3:4),3) + cytoFrac*u(i,:,5))./(yvals0 + cytoFrac*u(1,:,5)),'Color',cmap(floor(i/10)+1,:),'LineWidth',2); hold on;
    title('Fold-increase in signal (including cyto)');
end
xlim([0 max(x)]);
set(gca,'XTick',0:5:20);

%% Compute basic properties 
% Fast diffusion anchor
D=1;
iLIDconc=10.^[-2:0.05:1];
sspbConc=10.^[-2:0.05:1];
t=0:0.1:10;
t1=find(t==1);
t10=find(t==10);
inputFun=@(x) inputFunFromFit(x);

% structure for holding the simulation data
sim.t=t;
sim.pin0=p0;
sim.inputFun=inputFun;
sim.iLIDconc=iLIDconc;
sim.sspbConc=sspbConc;
sim.D=D;
sim.output=cell(length(iLIDconc),length(sspbConc));
sim.fracRec=nan(length(iLIDconc),length(sspbConc));
sim.maxRec=nan(length(iLIDconc),length(sspbConc));
sim.maxRecCytoLevel=nan(length(iLIDconc),length(sspbConc));
sim.basalRec=nan(length(iLIDconc),length(sspbConc));
sim.basalRecCytoLevel=nan(length(iLIDconc),length(sspbConc));
sim.xHalf1=nan(length(iLIDconc),length(sspbConc));
sim.xHalf10=nan(length(iLIDconc),length(sspbConc));
sim.tMax=nan(length(iLIDconc),length(sspbConc));

for i=1:length(iLIDconc)
    for j=1:length(sspbConc)
        pin=p0; 
        pin(3)=D;
        pin(4)=sspbConc(j);
        pin(5)=iLIDconc(i);
        [t,x,u]=iLID_PDE_model(t_range,pin,x_range,inputFun);
        output=sum(u(:,:,3:4),3);
        sim.output{i,j}=output;
        sim.fracRec(i,j)=max(mean(output.*x,2))/sum(x*sspbConc(j));  % Fraction of sspB recruited to the membrane
        [sim.maxRec(i,j),ind]=max(output(:,1));  %max recruitment at the cell center
        sim.maxRecCytoLevel(i,j)=u(ind,1,5);  %local cytoplasmic sspB concentration at time of max recruitment at the cell center
        sim.basalRec(i,j)=output(1,1);
        sim.basalRecCytoLevel(i,j)=u(1,1,5);
        sim.tMax(i,j)=mean(t(output(:,1)==max(output(:,1))));
        yvals=output - output(1,1);  %new recruitment, excluding baseline
        sim.xHalf1(i,j)=interp1(yvals(t1,:),x,0.5*max(yvals(t1,:)));
        sim.xHalf10(i,j)=interp1(yvals(t10,:),x,0.5*max(yvals(t10,:)));
    end
    fprintf('%i ',i);
end
sim1=sim;
fprintf('\n');
clear sim
%save('G:\Shared drives\Collins lab\general-matlab\modeling\iLID\PDE local recruitment modeling data 06-23-2020.mat','sim1');
%clear sim1

% **************************
%% Slow diffusing anchor
D=0.1;
iLIDconc=10.^[-2:0.05:1];
sspbConc=10.^[-2:0.05:1];
t=0:0.1:10;
t1=find(t==1);
t10=find(t==10);
inputFun=@(x) inputFunFromFit(x);
%inputFun0=@(x) gaussian(x,0,3.3) + 0.2*gaussian(x,0,6);
%inputFun=@(x) inputFun0(x)/inputFun0(0);
% structure for holding the simulation data
sim.t=t;
sim.pin0=p0;
sim.inputFun=inputFun;
sim.iLIDconc=iLIDconc;
sim.sspbConc=sspbConc;
sim.D=D;
sim.output=cell(length(iLIDconc),length(sspbConc));
sim.fracRec=nan(length(iLIDconc),length(sspbConc));
sim.maxRec=nan(length(iLIDconc),length(sspbConc));
sim.maxRecCytoLevel=nan(length(iLIDconc),length(sspbConc));
sim.basalRec=nan(length(iLIDconc),length(sspbConc));
sim.basalRecCytoLevel=nan(length(iLIDconc),length(sspbConc));
sim.xHalf1=nan(length(iLIDconc),length(sspbConc));
sim.xHalf10=nan(length(iLIDconc),length(sspbConc));
sim.tMax=nan(length(iLIDconc),length(sspbConc));

for i=1:length(iLIDconc)
    for j=1:length(sspbConc)
        pin=p0; 
        pin(3)=D;
        pin(4)=sspbConc(j);
        pin(5)=iLIDconc(i);
        [t,x,u]=iLID_PDE_model(t_range,pin,x_range,inputFun);
        output=sum(u(:,:,3:4),3);
        sim.output{i,j}=output;
        sim.fracRec(i,j)=max(mean(output.*x,2))/sum(x*sspbConc(j));  % Fraction of sspB recruited to the membrane
        [sim.maxRec(i,j),ind]=max(output(:,1));  %max recruitment at the cell center
        sim.maxRecCytoLevel(i,j)=u(ind,1,5);  %local cytoplasmic sspB concentration at time of max recruitment at the cell center
        sim.basalRec(i,j)=output(1,1);
        sim.basalRecCytoLevel(i,j)=u(1,1,5);
        sim.tMax(i,j)=mean(t(output(:,1)==max(output(:,1))));
        yvals=output - output(1,1);  %new recruitment, excluding baseline
        sim.xHalf1(i,j)=interp1(yvals(t1,:),x,0.5*max(yvals(t1,:)));
        sim.xHalf10(i,j)=interp1(yvals(t10,:),x,0.5*max(yvals(t10,:)));
    end
    fprintf('%i ',i);
end
sim01=sim;
fprintf('\nDone.\n');
clear sim;
%save('G:\Shared drives\Collins lab\general-matlab\modeling\iLID\PDE local recruitment modeling data 06-23-2020.mat','sim01','-append');
%clear sim01

%% save results
%save('G:\Shared drives\Collins lab\general-matlab\modeling\iLID\PDE local recruitment modeling data 06-18-2020.mat','sim01','sim1');
%% load results
%load('G:\Shared drives\Collins lab\general-matlab\modeling\iLID\PDE local recruitment modeling data 06-18-2020.mat','sim01','sim1');
load('G:\Shared drives\Collins lab\general-matlab\modeling\iLID\PDE local recruitment modeling data 06-23-2020.mat','sim01','sim1');

%% save results - slimmed down version, lacking the full simulation results
sim01=rmfield(sim01,'output');
sim1=rmfield(sim1,'output');
save('G:\Shared drives\Collins lab\general-matlab\modeling\iLID\PDE local recruitment modeling data slim 06-23-2020_kOffLit-0.05.mat','sim01','sim1');
%% load results - slimmed down version, lacking the full simulation results
load('G:\Shared drives\Collins lab\general-matlab\modeling\iLID\PDE local recruitment modeling data slim 06-23-2020.mat','sim01','sim1');

%% Figure for fast diffusion
cytoFrac=0.075;   % Fraction of cytosolic molecules visible in TIRF

sspbConc=sim1.sspbConc;
iLIDconc=sim1.iLIDconc;
tickInd=1:20:length(sspbConc);
figure;
subplot(2,4,1);
imagesc(sim1.tMax); colorbar;
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
xlabels=arrayfun(@(x) {num2str(x)},sspbConc(tickInd));
ylabels=arrayfun(@(x) {num2str(x)},iLIDconc(tickInd));
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Time to Max Recruitment')

subplot(2,4,2);
imagesc(sim1.xHalf1,[4 9]); 
colorbar('Ticks',4:8,'TickLabels',4:8);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Half-max spread at 1 sec');

subplot(2,4,3);
imagesc(sim1.xHalf10,[4 9]);
colorbar('Ticks',4:8,'TickLabels',4:8);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Half-max spread at 10 sec');

subplot(2,4,4);
imagesc(log2(sim1.maxRec-sim1.basalRec),log2([5e-4 5])); 
colorbar('Ticks',log2([1e-3 1e-2 0.1 1]),'TickLabels',[1e-3 1e-2 0.1 1]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Max Absolute Recruitment')

subplot(2,4,5);
imagesc(log2(sim1.basalRec),log2([5e-4 5])); 
colorbar('Ticks',log2([1e-3 1e-2 0.1 1]),'TickLabels',[1e-3 1e-2 0.1 1]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Basal Recruitment');

subplot(2,4,6);
imagesc(log2(sim1.maxRec./sim1.basalRec),log2([1 25]));
colorbar('Ticks',log2([1 2 5 10 20]),'TickLabels',[1 2 5 10 20]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Fold Enrichment by Recruitment');

subplot(2,4,7);
imagesc(log2((sim1.maxRec + cytoFrac*sim1.maxRecCytoLevel)./(sim1.basalRec + cytoFrac*sim1.basalRecCytoLevel)),log2([1 10]));
colorbar('Ticks',log2([1 2 5 10]),'TickLabels',[1 2 5 10]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Apparent Log2 Fold Enrichment with cytosolic contribution');
%% Figure for slow diffusion
cytFrac=0.075;   % Fraction of cytosolic molecules visible in TIRF

sspbConc=sim01.sspbConc;
iLIDconc=sim01.iLIDconc;
tickInd=1:20:length(sspbConc);
figure;
subplot(2,4,1);
imagesc(sim01.tMax); colorbar;
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
xlabels=arrayfun(@(x) {num2str(x)},sspbConc(tickInd));
ylabels=arrayfun(@(x) {num2str(x)},iLIDconc(tickInd));
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Time to Max Recruitment')

subplot(2,4,2);
imagesc(sim01.xHalf1,[4 9]); 
colorbar('Ticks',4:8,'TickLabels',4:8);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Half-max spread at 1 sec');

subplot(2,4,3);
imagesc(sim01.xHalf10,[4 9]);
colorbar('Ticks',4:8,'TickLabels',4:8);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Half-max spread at 10 sec');

subplot(2,4,4);
imagesc(log2(sim01.maxRec-sim01.basalRec),log2([5e-4 5])); 
colorbar('Ticks',log2([1e-3 1e-2 0.1 1]),'TickLabels',[1e-3 1e-2 0.1 1]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Max Absolute Recruitment')

subplot(2,4,5);
imagesc(log2(sim01.basalRec),log2([5e-4 5])); 
colorbar('Ticks',log2([1e-3 1e-2 0.1 1]),'TickLabels',[1e-3 1e-2 0.1 1]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Basal Recruitment');

subplot(2,4,6);
imagesc(log2(sim01.maxRec./sim01.basalRec),log2([1 25]));
colorbar('Ticks',log2([1 2 5 10 20]),'TickLabels',[1 2 5 10 20]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Fold Enrichment by Recruitment');

subplot(2,4,7);
imagesc(log2((sim01.maxRec + cytoFrac*sim01.maxRecCytoLevel)./(sim01.basalRec + cytoFrac*sim01.basalRecCytoLevel)),log2([1 10]));
colorbar('Ticks',log2([1 2 5 10]),'TickLabels',[1 2 5 10]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Apparent Log2 Fold Enrichment with cytosolic contribution');

%% Compare "optimal" regions in concentration space
figure;
mask1=sim1.maxRec./sim1.basalRec > 10 & sim1.maxRec > 0.1 & sim1.basalRec<0.1;
mask01=sim01.maxRec./sim01.basalRec > 10 & sim01.maxRec > 0.1 & sim01.basalRec<0.1;
imagesc(mask01*2 + mask1);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Optimal concentration regime');

%% Compare fold-recruitment depending on diffusion coefficients
figure;

subplot(1,2,1);
imagesc(log2(sim1.maxRec./sim1.basalRec),[0 4.5]); 
colorbar('Ticks',log2([1 2 5 10 20]),'TickLabels',[1 2 5 10 20]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Fold Recruitment (D=1)');

subplot(1,2,2);
imagesc(log2(sim01.maxRec./sim01.basalRec),[0 4.5]); 
colorbar('Ticks',log2([1 2 5 10 20]),'TickLabels',[1 2 5 10 20]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Fold Recruitment (D=0.1)');

%% Compare actual fold-recruitment with apparent fold-recruitment estimated from TIRF
%% load results from ODE model
load('G:\Shared drives\Collins lab\general-matlab\modeling\iLID\ODE global recruitment modeling data 06-19-2020.mat','sim');

%% Compare fold-recruitment depending on diffusion coefficients
figure;

subplot(1,3,1);
imagesc(log2(sim.absRec./sim.basalRec),[0 4.5]); 
colorbar('Ticks',log2([1 2 5 10 20]),'TickLabels',[1 2 5 10 20]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Fold Recruitment (ODE)');

subplot(1,3,2);
imagesc(log2(sim1.maxRec./sim1.basalRec),[0 4.5]); 
colorbar('Ticks',log2([1 2 5 10 20]),'TickLabels',[1 2 5 10 20]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Fold Recruitment (D=1)');

subplot(1,3,3);
imagesc(log2(sim01.maxRec./sim01.basalRec),[0 4.5]); 
colorbar('Ticks',log2([1 2 5 10 20]),'TickLabels',[1 2 5 10 20]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Fold Recruitment (D=0.1)');

%% Compare slow diffusion model to ODE model
imagesc(log2(sim01.maxRec./sim01.basalRec) - log2(sim.absRec./sim.basalRec))
colorbar('Ticks',log2([1 2 4]),'TickLabels',[1 2 4]);
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Log2 Fold Increase in Recruitment (D=0.1 vs ODE)');

%% Make contour plot overlaying two different output parameters
cytoFrac=0.075;
figure('Position',[300 400 1200 500]);
subplot(1,2,1);
contour(sim1.maxRec./sim1.basalRec,[2 5 10 20],'r'); hold on;
contour(sim1.maxRec./sim1.basalRec,[10 10],'r','LineWidth',2); hold on;  %,'ShowText','on'
%contour((sim1.maxRec + cytoFrac*sim1.maxRecCytoLevel)./(sim1.basalRec + cytoFrac*sim1.basalRecCytoLevel),[2 5 10 20],'g','ShowText','on');
hold on;
contour(sim1.maxRec,[0.05 0.1 0.2 0.5 1],'LineColor',[0.85 0.85 0.85]);
contour(sim1.maxRec,[0.1 0.1],'k','LineWidth',2);
contour(sim1.basalRec,[0.01 0.02 0.05 0.1 0.2],'LineColor',[0.8 0.8 1]);
contour(sim1.basalRec,[0.05 0.05],'b','LineWidth',2);
axis square
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Contour Plot D=1');

subplot(1,2,2);
contour(sim01.maxRec./sim01.basalRec,[2 5 10 20],'r'); hold on;
contour(sim01.maxRec./sim01.basalRec,[10 10],'r','LineWidth',2); hold on;  %,'ShowText','on'
%contour((sim01.maxRec + cytoFrac*sim01.maxRecCytoLevel)./(sim01.basalRec + cytoFrac*sim01.basalRecCytoLevel),[2 5 10 20],'g','ShowText','on');
hold on;
contour(sim01.maxRec,[0.05 0.1 0.2 0.5 1],'LineColor',[0.85 0.85 0.85]);
contour(sim01.maxRec,[0.1 0.1],'k','LineWidth',2);
contour(sim01.basalRec,[0.01 0.02 0.05 0.1 0.2],'LineColor',[0.8 0.8 1]);
contour(sim01.basalRec,[0.05 0.05],'b','LineWidth',2);
axis square
xlabel('sspB Conc');
ylabel('iLID Conc');
set(gca,'XTick',tickInd);
set(gca,'YTick',tickInd);
set(gca,'YDir','normal');
set(gca,'XTickLabel',xlabels);
set(gca,'YTickLabel',ylabels);
title('Contour Plot D=0.1');

%% Plot to compare model to experimental data for local recruitment
cytoFrac=0.1;
figure;
obsFoldEnrich1=(sim1.maxRec + cytoFrac*sim1.maxRecCytoLevel)./(sim1.basalRec + cytoFrac*sim1.basalRecCytoLevel);
obsFoldEnrich01=(sim01.maxRec + cytoFrac*sim01.maxRecCytoLevel)./(sim01.basalRec + cytoFrac*sim01.basalRecCytoLevel);
ind=35; fprintf('sspB conc = %.2f\n',sspbConc(ind));
ind2=find(iLIDconc>=0.01 & iLIDconc<=2);
%scatter(sim1.tMax(ind2,ind),obsFoldEnrich(ind2,ind));
subplot(1,2,1);
plot(iLIDconc(ind2),obsFoldEnrich1(ind2,ind));
hold on;
plot(iLIDconc(ind2),obsFoldEnrich01(ind2,ind));
title('Observed Fold-Increase')
xlabel('iLID Conc (uM)');
legend({'D=1','D=0.1'});
ylim([1 6]);
subplot(1,2,2);
plot(iLIDconc(ind2),sim1.tMax(ind2,ind));
title('Time to Max Recruitment')
xlabel('iLID Conc (uM)');