close all;
clear; clc;

domain.xa   = 0;
domain.xb   = 2*pi;
domain.ya   = 0;
domain.yb   = 2*pi;
domain.za   = 0;
domain.zb   = 2*pi;

Lx = domain.xb - domain.xa;
Ly = domain.yb - domain.ya;
Lz = domain.zb - domain.za;

N  = 64;
Nx = N;
Ny = N;
Nz = N;

Plot = 2;

if Plot == 1 

%     scheme1 = 'linear';
    scheme1 = 'nonlinear';

%     scheme2 = '_1st';   % First-order scheme
%     scheme2 = '_2cn';   % Second-order CN scheme
    scheme2 = '_bdf2';  % Second-order BDF scheme

scheme = [scheme1 , scheme2];

else
    sch1 = 'linear';
%     sch1 = 'nonlinear';

    sch2 = '_1st'; 
%     sch2 = '_bdf2';

    sch3 = 'linear';
%     sch3 = 'nonlinear';

%     sch4 = '_1st'; 
    sch4 = '_bdf2';

    scheme1 = [sch1 , sch2];
    scheme2 = [sch3 , sch4];
end

% Parameters
para.epsilon = 6*pi/128;
para.M = 1; 
para.S1 = 4;
para.S2 = 4;
para.S3 = 1; 

PDE = 'data1';

if 1 == strcmp(PDE,'data1')
    T = 0.002;
    dt_array = 0.0001./2.^(1:6)';
    dt_ref = 1e-7;
    para.C0 = 100; 
    pde = ex03_1_Vesicles_data(para);
end

if 1 == strcmp(PDE,'data2')
    T = 0.002;
    dt_array = 0.001./2.^(6:11)';
    dt_ref = 1e-7;
    para.C0 = 10000;
    pde = ex03_2_Vesicles_data(para);
end

% Time: dt T

t0 = 0;
tsave = 0.2*T;

maxIt = length(dt_array);

if Plot == 1
error=zeros(maxIt,1);
order=zeros(maxIt,1);
if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
    name=[scheme,'_phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'dt=',num2str(dt_ref)];
    filename=[name '.mat'];
    load(filename,'phi');
    phi_exact = phi;
    clear phi;
else
    Lx = domain.xb - domain.xa;
    Ly = domain.yb - domain.ya;
    Lz = domain.zb - domain.za;
    hx = Lx/Nx;
    hy = Ly/Ny;
    hz = Lz/Nz;
    x  = domain.xa + hx*(0:Nx-1);
    y  = domain.ya + hy*(0:Ny-1);
    z  = domain.za + hz*(0:Nz-1);
    [xx,yy,zz] = ndgrid(x,y,z);
    phi_exact = pde.exact(xx,yy,zz,T);
end
for k = 1:maxIt
    dt = dt_array(k);
    name=[scheme,'_phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
        'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'dt=',num2str(dt)];
    filenamek=[name '.mat'];
    load(filenamek,'phi','hx','hy','hz');
    err    = fftn((phi_exact - phi).^2);
    error(k,1) = sqrt(err(1,1)*hx*hy*hz);   % L2
    clear phi;
end
order(2:maxIt) = log(error(1:maxIt-1)./error(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));

%% Plot(1)
figure(2)
hh=loglog(dt_array,error,'*-','LineWidth',3,'MarkerSize',10);
hold on 
tf = strcmp(scheme2,'_1st');
if tf == 1 
    hh=loglog(dt_array,dt_array,'k-.','LineWidth',3,'MarkerSize',10);
    hold on
    h = legend('First scheme: $\phi$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
else 
    hh=loglog(dt_array,dt_array.^2,'k-.','LineWidth',3,'MarkerSize',10);
    hold on
    h = legend('Second scheme: $\phi$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
end
set(h,'interpreter','latex','FontSize',15);
xlabel('Time step $\delta t$','Interpreter','latex');
%  ylabel('Error\_max')
ylabel('$L^2$ error','Interpreter','latex');
% ylim([1e-12 5e-3])
ylim([1e-14 5e-5])
grid on;
hold on;

figname1 = ['C:\Users\heyan\Desktop\some files\figs\accuracy\',scheme,'_accuracy_3d','.png'];
% print(figname1,'-dpng', '-r300')

%% Save error and order(1)
name=[scheme,'_phi_e',num2str(para.epsilon),'M',num2str(para.M),...
      'Nx=',num2str(N),'Ny=',num2str(N),'dt1=',num2str(dt)];
fileID = fopen([name,'.txt'],'w');
% fprintf(fileID,'%6s\n','%% Results');
% fprintf(fileID,'%6s\n','% dt	   &   Error_L2	   &  Order');
% A = [dt_array error];
% fprintf(fileID,'%.12f   %.4e   \n',A');
fprintf(fileID,'%.12f     %.4e      %.2f \n',[dt_array,error,order]');
fclose(fileID);

else
%% Compute order of convergence(2)
error1=zeros(maxIt,1);
order1=zeros(maxIt,1);
error2=zeros(maxIt,1);
order2=zeros(maxIt,1);
if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
    name1=[scheme1,'_phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'dt=',num2str(dt_ref)];
    name2=[scheme2,'_phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'dt=',num2str(dt_ref)];
    filename1=[name1 '.mat'];
    filename2=[name2 '.mat'];
    load(filename1,'phi');
    phi_exact1 = phi;
    clear phi
    load(filename2,'phi');
    phi_exact2 = phi;
    clear phi;
else
    Lx = domain.xb - domain.xa;
    Ly = domain.yb - domain.ya;
    Lz = domain.zb - domain.za;
    hx = Lx/Nx;
    hy = Ly/Ny;
    hz = Lz/Nz;
    x  = domain.xa + hx*(0:Nx-1);
    y  = domain.ya + hy*(0:Ny-1);
    z  = domain.za + hz*(0:Nz-1);
    [xx,yy,zz] = ndgrid(x,y,z);
    phi_exact = pde.exact(xx,yy,zz,T);
end
for k = 1:maxIt
    dt = dt_array(k);
    name1=[scheme1,'_phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'dt=',num2str(dt)];
    name2=[scheme2,'_phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'dt=',num2str(dt)];
    filenamek1=[name1 '.mat'];
    load(filenamek1,'phi','hx','hy','hz');
    err1    = fftn((phi_exact1 - phi).^2);
    error1(k,1) = sqrt(err1(1,1,1)*hx*hy*hz);   % L2
    clear phi;
    filenamek2=[name2 '.mat'];
    load(filenamek2,'phi','hx','hy','hz');
    err2    = fftn((phi_exact2 - phi).^2);
    error2(k,1) = sqrt(err2(1,1,1)*hx*hy*hz);   % L2
    clear phi;
end
order1(2:maxIt) = log(error1(1:maxIt-1)./error1(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));
order2(2:maxIt) = log(error2(1:maxIt-1)./error2(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));


%% Plot(2)
figure(2)
hh=loglog(dt_array,error1,'*-','LineWidth',2,'MarkerSize',10);
hold on
hh=loglog(dt_array,error2,'p-','LineWidth',2,'MarkerSize',10);
hold on 
tf1 = strcmp(sch2,'_1st');
tf2 = strcmp(sch4,'_1st');
if tf1 == 1 && tf2 == 1
    hh=loglog(dt_array,dt_array,'k-.','LineWidth',3,'MarkerSize',10);
    hold on
    h = legend('First scheme: $\phi$','Second scheme: $\phi$','$\mathcal{O}(\delta t)$','Location','SouthEast');
elseif tf1 == 0 && tf2 == 0
    hh=loglog(dt_array,dt_array.^2,'k-.','LineWidth',3,'MarkerSize',10);
    hold on 
    h = legend('First scheme: $\phi$','Second scheme: $\phi$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
else 
    hh=loglog(dt_array,dt_array,'k--','LineWidth',3,'MarkerSize',10);
    hold on
    hh=loglog(dt_array,dt_array.^2,'k-.','LineWidth',3,'MarkerSize',10);
    hold on
    h = legend(' First scheme: $\phi$','First scheme: $\phi$','$\mathcal{O}(\delta t)$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
end
set(h,'interpreter','latex','FontSize',15);
xlabel('Time step $\delta t$','Interpreter','latex');
%  ylabel('Error\_max')
ylabel('$L^2$ error','Interpreter','latex');
ylim([1e-14 1e-4])
% ylim([1e-7 5e-4])
grid on;
hold on;

set(gca,'FontSize',18);
set(gca,'linewidth',1.1)

scheme = [scheme1, scheme2];

figname1 = ['C:\Users\heyan\Desktop\some files\figs\accuracy\',scheme1,'and',scheme2,'_accuracy_3d','.png'];
print(figname1,'-dpng', '-r300')

end