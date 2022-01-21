
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

%% Plot experimental timecourses
% Rearrange data for center point in a matrix of cell x time
for i=1:4
    datMat{i}=nan(size(yOutFinal{i}{1},1),length(yOutFinal{i}));
    for j=1:length(yOutFinal{i})
        datMat{i}(:,j)=yOutFinal{i}{j}(:,1);
    end
end

%% Plot timecourse for center region
kRevert = 0.02;  %iLID inactivation rate
kOffLit = 0.2;  %iLID-SspB association rate in the dark state
kBind = kOffLit/kDlight *2;  %iLID-SspB association rate
kOffDark = kBind*kDdark;  %iLID-SspB association rate in the active/lit state
SspBTot = 0.25;  %total concentration of SspB
iLIDTot = 0.5;   %total concentration of iLID
D=1;   %iLID membrane lateral diffusion rate

t_range=0:0.1:20;

p0=[kRevert kBind D SspBTot iLIDTot cellRadius];
[t,x,u]=iLID_PDE_model(t_range,p0,x_range,inputFunFromFit);
plot(t,sum(u(:,1,3:4),3));
