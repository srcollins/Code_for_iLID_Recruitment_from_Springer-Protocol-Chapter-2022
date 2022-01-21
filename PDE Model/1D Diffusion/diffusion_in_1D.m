function [t,x,u]=diffusion_in_1D(t_range,x_range,p,u0)
% this function simulates 1-dimensional diffusion along a line (such as a
% cell edge)

%% Define constants

global D;
global domain_length;

D=p(1);    %in units of microns squared per second, based on estimates from Jim Ferrell
domain_length=p(2);   %cell radius in microns
n_pts=length(u0);

%% Set up region, bounds

m = 0;      %symmetry, 0=slab, 1=cylindrical, 2=spherical
step=(x_range(end)-x_range(1))/(n_pts-1);
x = [(-domain_length/2):step:(x_range(1)-step) x_range (x_range(end)+step):step:(domain_length/2)];  
t = t_range;
temp=zeros(size(x));
temp(x<x_range(1))=1;
temp(x>x_range(end))=1;
temp(x>=x_range(1) & x<=x_range(end))=u0;
u0=temp;

%% Do the computation

%sol = pdepe(m,@diff1Dpde,@diff1Dic,@diff1Dbc,x,t);
size(x);
size(u0);
sol = pdepe(m,@diff1Dpde,@(z) interp1(x,u0,z),@diff1Dbc,x,t);
% Extract the first solution component as u.  This is not necessary
% for a single equation, but makes a point about the form of the output.
u = sol(:,:,1);

%% Display results

% A surface plot is often a good way to study a solution.
figure;
surf(x,t,u);    
title('Numerical solution');
xlabel('Distance x');
ylabel('Time t');

% A solution profile can also be illuminating.
figure;
plot(x,u(end,:),'o');
title(['Solutions at t = ' num2str(t(end)) '.']);
%legend('Numerical, 20 mesh points',0);
xlabel('Distance x');
ylabel('u(x,2)');

% --------------------------------------------------------------------------
function [c,f,s] = diff1Dpde(x,t,u,DuDx)
global D;

c = 1;
f = D * DuDx;
s = 0;

% --------------------------------------------------------------------------
function u0 = diff1Dic(x)
global cell_radius;

u0 = double(x==cell_radius);
    
% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = diff1Dbc(xl,ul,xr,ur,t)
pl = 0;
ql = 1;
pr = 0;
qr = 1;
