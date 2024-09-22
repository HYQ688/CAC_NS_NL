clear;
clc;
clf;
close all;

dirname = 'ex20_BDF2_3D_CAC_Vesicles_data_LM_four_globe';


datadir = [dirname,'/data'];
figdir  = [dirname,'/'];

% Space: Domain and N
domain.xa =   0;
domain.xb = 2*pi;
domain.ya =   0;
domain.yb = 2*pi;
domain.za =   0;
domain.zb = 2*pi;

% N = 96;
N = 64;
Nx = N; Ny = N; Nz = N; 
Lx = domain.xb - domain.xa;
Ly = domain.yb - domain.ya;
Lz = domain.zb - domain.za;
hx = Lx/Nx;
hy = Ly/Ny;
hz = Ly/Nz;
x  = domain.xa + hx*(0:Nx-1);
y  = domain.ya + hy*(0:Ny-1);
z  = domain.za + hy*(0:Nz-1);

[xx,yy,zz] = ndgrid(x,y,z);

tsave = 0.1;
kk = 0;
% figure(4) 
% for t= [0:0.04:0.76 1.44 1.6 1.8 2]
for t = 0:0.04:2
% for t = 0
t 

      %% figures u
%     kk = kk+1;
%     figure(kk)
figure(1)
    clf
%     filename = [datadir '/phi_t=' num2str(t)];
%     ss = [filename '.txt'];
%     phiu = load(ss);
%     filename = [datadir '/phi_b_t=' num2str(t)];
%     tt = [filename '.txt'];
%     phiv = load(tt);
%     phiu = reshape(phiu-phiv,Nx,Ny,Nz);    
%     p0 = patch(isosurface(xx,yy,zz,phiu,0.001));
%     set(p0,'FaceColor','red', 'EdgeColor','none');
%     hold on
%     p0 = patch(isosurface(xx,yy,zz,phiu,-0.001));
%     set(p0,'FaceColor','blue', 'EdgeColor','none');
    

    filename = [datadir '/phi_t=' num2str(t)];
    ss = [filename '.txt'];
    phiu = load(ss);
    phiu = reshape(phiu,Nx,Ny,Nz);    
    p0 = patch(isosurface(xx,yy,zz,phiu,0));
    set(p0,'FaceColor','red', 'EdgeColor','none');

%     hold on
% 
%     filename = [datadir '/phi_b_t=' num2str(t)];
%     ss = [filename '.txt'];
%     phiu = load(ss);
%     phiu = reshape(phiu,Nx,Ny,Nz);    
%     p0 = patch(isosurface(xx,yy,zz,phiu,0.2));
%     set(p0,'FaceColor','blue', 'EdgeColor','none');

    daspect([1 1 1]); 
    camlight;  
    lighting phong;
%     axis image;
%     axis tight;
    axis square;
    axis equal;
%     axis([min(xx(:)) min(xx(:))min(yy(:)) max(yy(:)) min(zz(:)) max(zz(:)) ]);
%     axis 'auto xy';
    xlim([0, 2*pi])
    ylim([0, 2*pi])
    zlim([0, 2*pi])
    box on; 

    view(-161,16);
%     view(115,20)    
    
    ax = gca;
    c = ax.Color;
    ax.Color = [0.30,0.75,0.90];

%     xlabel('$x$','FontSize',18,'Interpreter','LaTex');
%     ylabel('$y$','FontSize',18,'Interpreter','LaTex');
%     zlabel('$z$','FontSize',18,'Interpreter','LaTex');
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     set(gca,'ztick',[])
    drawnow;
    set(gcf, 'InvertHardCopy', 'off');
    set(0,'defaultfigurecolor','w') 
    figname = ['C:\Users\Administrator\Desktop\论文\figs\3D\18a\',dirname,'_phi_t=' num2str(t) '.png'];
%     print(figname,'-dpng', '-r300')
   
end
