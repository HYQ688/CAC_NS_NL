close all;
clear; clc;

% pathname1 = 'D:\学校文件\论文\code\00\New_LM_p_SAV\CAC_Vesicle_newLM_p_SAV\CAC_Vesicles_surface_LagrangeMultiplier_SAV_ok\example02_2D_First_order_linear_scheme_in_Time_reference_ok';
% pathname2 = 'D:\学校文件\论文\code\00\New_LM_p_SAV\CAC_Vesicle_newLM_p_SAV\CAC_Vesicles_surface_LagrangeMultiplier_SAV_ok\example03_2D_First_order_Order_in_Time_reference_ok';
pathname1 = 'D:\学校文件\论文\code\00\New_LM_p_SAV\CAC_Vesicle_newLM_p_SAV\CAC_Vesicles_surface_LagrangeMultiplier_SAV_ok\example08_2D_BDF2_linear_scheme_in_Time_reference_ok';
pathname2 = 'D:\学校文件\论文\code\00\New_LM_p_SAV\CAC_Vesicle_newLM_p_SAV\CAC_Vesicles_surface_LagrangeMultiplier_SAV_ok\example09_2D_BDF2_Order_in_Time_reference_ok';




domain.left   = 0;
domain.right  = 2*pi;
domain.bottom = 0;
domain.top    = 2*pi;

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;

N  = 64;
Nx = N;
Ny = N;

% Parameters
para.C0 = 10000; % SAV
para.epsilon = 6*pi/128;
para.M = 1;

addpath( genpath(pathname1) );   

PDE = 'data1';

if 1 == strcmp(PDE,'data1')
    T = 0.0002;
    dt_array = 0.00001./2.^(0:6)';
    dt_ref = 1e-8;
    pde = ex03_1_Vesicles_data(para);
end

if 1 == strcmp(PDE,'data2')
    T = 0.0001;
    dt_array = 0.0001./2.^(5:16)';
    dt_ref = 1e-9;
    pde = ex03_2_Vesicles_data(para);
end

% Time: dt T

t0 = 0;
tsave = 0.2*T;

maxIt = length(dt_array);

%% Compute order of convergence
error1=zeros(maxIt,1);
order1=zeros(maxIt,1);
if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
    name=['phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt_ref)];
    filename=[name '.mat'];
    load(filename,'phi');
    phi_exact = phi;
    clear phi;
else
    Lx = domain.right - domain.left;
    Ly = domain.top   - domain.bottom;
    hx = Lx/Nx;
    hy = Ly/Ny;
    x  = domain.left   + hx*(0:Nx-1);
    y  = domain.bottom + hy*(0:Ny-1);
    [xx,yy] = ndgrid(x,y);
    phi_exact = pde.exact(xx,yy,T);
end
for k = 1:maxIt
    dt = dt_array(k);
    name=['phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
        'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];
    filenamek=[name '.mat'];
    load(filenamek,'phi','hx','hy');
    err    = fft2((phi_exact - phi).^2);
    error1(k,1) = sqrt(err(1,1)*hx*hy);   % L2
    clear phi;
end
order1(2:maxIt) = log(error1(1:maxIt-1)./error1(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));

rmpath( genpath(pathname1)  );

addpath( genpath(pathname2) );  

PDE = 'data1';

if 1 == strcmp(PDE,'data1')
    T = 0.0002;
    dt_array = 0.00001./2.^(0:6)';
    dt_ref = 1e-8;
    pde = ex03_1_Vesicles_data(para);
end

if 1 == strcmp(PDE,'data2')
    T = 0.0001;
    dt_array = 0.0001./2.^(5:16)';
    dt_ref = 1e-9;
    pde = ex03_2_Vesicles_data(para);
end

% Time: dt T

t0 = 0;
tsave = 0.2*T;

maxIt = length(dt_array);

%% Compute order of convergence
error2=zeros(maxIt,1);
order2=zeros(maxIt,1);
if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
    name=['phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt_ref)];
    filename=[name '.mat'];
    load(filename,'phi');
    phi_exact = phi;
    clear phi;
else
    Lx = domain.right - domain.left;
    Ly = domain.top   - domain.bottom;
    hx = Lx/Nx;
    hy = Ly/Ny;
    x  = domain.left   + hx*(0:Nx-1);
    y  = domain.bottom + hy*(0:Ny-1);
    [xx,yy] = ndgrid(x,y);
    phi_exact = pde.exact(xx,yy,T);
end
for k = 1:maxIt
    dt = dt_array(k);
    name=['phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
        'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];
    filenamek=[name '.mat'];
    load(filenamek,'phi','hx','hy');
    err    = fft2((phi_exact - phi).^2);
    error2(k,1) = sqrt(err(1,1)*hx*hy);   % L2
    clear phi;
end
order2(2:maxIt) = log(error2(1:maxIt-1)./error2(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));

rmpath( genpath(pathname2)  );


%% Plot
figure(2)
hh=loglog(dt_array,error1,'*-','LineWidth',2,'MarkerSize',10);
hold on
hh=loglog(dt_array,error2,'p-','LineWidth',2,'MarkerSize',10);
hh=loglog(dt_array,dt_array.^2,'k-.','LineWidth',3,'MarkerSize',10);
h = legend('Linear scheme: $\phi$','Non linear scheme: $\phi$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
set(h,'interpreter','latex','FontSize',15);
xlabel('Time step $\delta t$','Interpreter','latex');
%  ylabel('Error\_max')
ylabel('$L^2$ error','Interpreter','latex');
grid on;
hold on;

set(gca,'FontSize',22);
set(gca,'linewidth',1.1)

% figname1 = ['C:\Users\heyan\Desktop\sjj\songfig\error_first_order','.png'];
% print(figname1,'-dpng', '-r300')

figname1 = ['C:\Users\heyan\Desktop\sjj\songfig\error_BDF2','.png'];
print(figname1,'-dpng', '-r300')