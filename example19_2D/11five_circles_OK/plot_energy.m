clear;
clc;
clf;
close all;

%     sch1 = 'linear';
%     sch1 = 'nonlinear';
    sch1 = 'MSAV';

%     sch2 = '_1st'; 
    sch2 = '_bdf2';

    sch3 = 'linear';
%     sch3 = 'nonlinear';
%     sch3 = 'MSAV';

%     sch4 = '_1st'; 
    sch4 = '_bdf2';

    sch5 = 'nonlinear';
    sch6 = '_bdf2';
    scheme1 = [sch1 , sch2];
    scheme2 = [sch3 , sch4];
    scheme3 = [sch5,sch6];


pdename1 = [scheme1,'_ex03_Vesicles_data';];
pdename2 = [scheme2,'_ex02_Vesicles_data';];
pdename3 = [scheme3,'_ex02_Vesicles_data';];



N  = 64;
Nx = N; 
Ny = N;
domain.left   =  0;
domain.right  = 2*pi;
domain.bottom =  0;
domain.top    = 2*pi;
Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;
hx = Lx/Nx;
hy = Ly/Ny;

x  = domain.left   + hx*(0:Nx-1);
y  = domain.bottom + hy*(0:Ny-1);

[xx,yy] = ndgrid(x,y);

[k_x,k_y,kx,ky,kxx,kyy,k2,k4] = prepare_fft2_v2(Lx,Ly,Nx,Ny);

dt_ref = 0.01./(2.^4);

hold on

lineType ={'k--','b:','r-.','g:.'};

% figure(1)
% % subplot(2,1,1)
% %% energy 
% energy=load([pdename,num2str(dt_ref),'_energy.txt']);
% tmp = 1;
% yyaxis left; % 激活左边的轴
% plot(energy(tmp:1:end,1),energy(tmp:1:end,3),'-','LineWidth',4.5);
% xlabel('Time','Fontsize',20);ylabel('Energy $E_{bdf}$','Fontsize',20,'interpreter','latex');
% xlim([0,4])
% % ylim([-0.1,1.4])
% % % yticks([0:2:20])
% set(gca,'FontSize',22);
% set(gca,'linewidth',1.8)
% % set(gca,'xtick',0:0.2:T);
% % set(gca,'ytick',3:0.1:4);
% % grid on;
% box on;
% figure_FontSize=35;
% % set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% % set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(findobj('FontSize',10),'FontSize',figure_FontSize); %这4句是将字体大小改为8号字，在小图里很清晰
% 
% figname1 = ['../../../papers/paper15_Vesicles/01_CHNS_cn/figure_Vesicles/',dirname,'_stability','.png'];
% % print(figname1,'-dpng', '-r300')
% 
% 
% 
% % figure(2)
% % subplot(2,1,2)
%% mass
mass1=load([pdename1,num2str(dt_ref),'_mass.txt']);
mass2=load([pdename2,num2str(dt_ref),'_mass.txt']);
mass3=load([pdename3,num2str(dt_ref),'_mass.txt']);
tmp = 1;
% yyaxis right; % 激活右边的轴
   fillcolor1=[0.85, 0.33, 0.10];
   fillcolor2=[0.93, 0.69, 0.13];
   fillcolor3=[0.00, 0.45, 0.74];
plot(mass1(tmp:1:end,1),abs((mass1(tmp:1:end,3)-mass1(1,3))./mass1(1,3)),'-','LineWidth',2.5,'Color',fillcolor3);
hold on
plot(mass2(tmp:1:end,1),abs((mass2(tmp:1:end,3)-mass2(1,3))./mass2(1,3)),'-','LineWidth',2.5,'Color',fillcolor1);
hold on
plot(mass3(tmp:1:end,1),abs((mass3(tmp:1:end,3)-mass3(1,3))./mass3(1,3)),'-','LineWidth',2.5,'Color',fillcolor2);
xlabel('Time','Fontsize',20);ylabel('Surface','Fontsize',20,'interpreter','latex');
xlim([0,2])
% if 1 == strcmp(dirname,'ex17_Vesicles_data_M2_50')
%     ylim([16,16.8])
% end
% % yticks([0:2:20])
% set(gca,'FontSize',22);
% set(gca,'linewidth',1.8)
% set(gca,'xtick',0:0.2:T);
% set(gca,'ytick',3:0.1:4);
% grid on;
box on;
% figure_FontSize=35;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize);
% set(get(gca,'YLabel'),'FontSize',figure_FontSize);
% set(findobj('FontSize',10),'FontSize',figure_FontSize); %这4句是将字体大小改为8号字，在小图里很清晰
legend('MSAV','NLM0','NLM3','Fontsize',20,'interpreter','latex','box','off');

figname2 = ['C:\Users\heyan\Desktop\some files\figs\deformation\2d50\',scheme1,'_and_',scheme2,'_and_',scheme3,'_stability_surface_2D','.png'];
print(figname2,'-dpng', '-r300')