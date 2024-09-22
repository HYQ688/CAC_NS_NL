function [phi, u0, v0] = CAC_Vesicle_2D_MSAV_1st(pde,domain,Nx,Ny,time,option)
% Solve 2D Vesicle model
% Qi Li
% Reference: Li, Xi, Tongmao Li, Rungting Tu, Kejia Pan, Chuanjun Chen, and Xiaofeng Yang. 
% "Efficient energy stable scheme for volume-conserved phase-field elastic bending energy model of lipid vesicles." 
% Journal of Computational and Applied Mathematics 385 (2021): 113177.
% 04/09/2022
global dt kx ky kxx kyy k2 k4 hx hy Lx Ly ...
       epsilon M beta beta_m C0 M1 M2 S1 S2 S3

if ~exist('option','var'), option = []; end
if ~isfield(option,'tol')
    option.tol = 10^-14;   % default tol
end
if ~isfield(option,'tolit')
    option.tolit = 10^-8;   % default tolit
end
if ~isfield(option,'maxit')
    option.maxit = 2000;   % default maxit
end
if ~isfield(option,'plotflag')
    option.plotflag = 0;   
end
if ~isfield(option,'saveflag')
    option.saveflag = 0;  
end
if ~isfield(option,'savefinal')
    option.savefinal = 0;  
end
if ~isfield(option,'printflag')
    option.printflag = 0;   
end
if ~isfield(option,'vtkflag')
    option.printflag = 0;   
end
if ~isfield(option,'energyflag')
    option.energyflag = 0;   
end
if 1 == option.energyflag
    figname_mass = [pde.name,num2str(time.dt),'_mass.txt'];
    figname_energy = [pde.name,num2str(time.dt),'_energy.txt'];      
    out1 = fopen(figname_mass,'w');
    out2 = fopen(figname_energy,'w');
end

tol = option.tol;
tolit = option.tolit;
maxit = option.maxit;

%%
T  = time.T;
t  = time.t0;
dt = time.dt;
tsave = time.tsave;

dir_fig  = [pde.name '/fig'];
dir_data = [pde.name '/data'];

epsilon = pde.epsilon;
M       = pde.M;
C0      = pde.C0;
beta_m  = pde.beta_m;
M1   = pde.M1;
M2   = pde.M2;
S1   = pde.S1;
S2   = pde.S2;
S3   = pde.S3;

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;
hx = Lx/Nx;
hy = Ly/Ny;
% x  = domain.left   + hx*(1:Nx);
% y  = domain.bottom + hy*(1:Ny);
x  = domain.left   + hx*(0:Nx-1);
y  = domain.bottom + hy*(0:Ny-1);

% [k_x,k_y,kx,ky,kxx,kyy,k2,k4] = prepare_fft2_v2(Lx,Ly,Nx,Ny);
k_x = 1i*[0:Nx/2 -Nx/2+1:-1]*(2*pi/Lx);
k_y = 1i*[0:Ny/2 -Ny/2+1:-1]*(2*pi/Ly);
[kx, ky] = ndgrid(k_x,k_y);

k2x = k_x.^2;
k2y = k_y.^2;
[kxx, kyy] = ndgrid(k2x,k2y);
k2 = kxx + kyy;
k4 = k2.^2;

[xx,yy] = ndgrid(x,y);
phi0 = pde.init(xx,yy);
nfigure =1;

beta  = beta_m*B(phi0);

%% plot initial value
nfigure =1;
if 1 == option.saveflag
    if ~exist(dir_data,'dir')
        mkdir(dir_data);
    end
    ss1 = [dir_data '/phi_t=' num2str(t) '.txt'];
    writematrix(phi0,ss1,'Delimiter',' ');
    writematrix(xx,[dir_data '/X.txt'],'Delimiter',' ');
    writematrix(yy,[dir_data '/Y.txt'],'Delimiter',' ');
end
if 1 == option.plotflag
    if 1 == option.saveflag
        showsolution_2D(nfigure,xx,yy,phi0,t,dir_fig);
    else
        showsolution_2D(nfigure,xx,yy,phi0,t);
    end
end


nplot = round((T-t)/dt);
nsave = round(tsave/dt);

tstart = tic;

u0 = fun_u_init(phi0);
v0 = fun_v_init(phi0);

% Initial energy
if 1 == option.energyflag
    calculate_energy(out1,out2,t,phi0,u0,v0);
end

for nt = 1:nplot
    t = t+dt;
     
    phi_star = phi0;    
    
    % step 1
    H = fun_H(phi_star);
    G = fun_G(phi_star);
    g = get_rhs(phi0,u0,v0,H,G);
    if isfield(pde,'rhs') && isfield(pde,'exact')
        ephi   = pde.exact(xx,yy,t);
        ephi_t = pde.exact_t(xx,yy,t);
        B_phi = fft2(epsilon/2.*grad_square(ephi)+1.0./epsilon.*F(ephi));
        B_phi = B_phi(1,1).*hx*hy;
        emu   = -epsilon.*lap_diff(omega(ephi)) + 1/epsilon.*f_der(ephi).*omega(ephi)...
                +M2.*(B_phi-beta).*epsilon.*omega(ephi);
        emu_bar = fft2(emu);
        emu_bar = emu - emu_bar(1,1) / (Lx*Ly);
        tmp  = (ephi_t./M + emu_bar) / epsilon;

%         tmp = pde.rhs(xx,yy,t,beta);
%         subplot(1,2,1)
%         mesh(tmp)
%         subplot(1,2,2)
%         mesh(tmp2)
        g = g + tmp;
    end

    psiA = inv_A(H);
    psiB = inv_A(G);
    psiC = inv_A(g);  
    
    b1 = fft2(H.*psiC);
    b1 = b1(1,1)*hx*hy;
    
    b2 = fft2(G.*psiC);
    b2 = b2(1,1)*hx*hy; 
    
    mu1 = fft2(H.*psiA);
    mu1 = mu1(1,1)*hx*hy;
    
    mu2 = epsilon.*M2*fft2(H.*psiB);
    mu2 = mu2(1,1)*hx*hy;
    
    mu3 = fft2(G.*psiA);
    mu3 = mu3(1,1)*hx*hy;
    
    mu4 = epsilon.*M2*fft2(G.*psiB);
    mu4 = mu4(1,1)*hx*hy;
    
    % Step 2
    X = (b1 - b2*mu2 + b1*mu4)./(mu1 + mu4 + mu1*mu4 - mu2*mu3 + 1);
    Y = (b2 + b2*mu1 - b1*mu3)./(mu1 + mu4 + mu1*mu4 - mu2*mu3 + 1);
   
    % Step 3
    phi = -X.*psiA - epsilon.*M2*Y.*psiB + psiC; 
    
    %% update phi0
    u0 = fun_u(phi,phi0,u0,H); 
    v0 = fun_v(phi,phi0,v0,G); 
    phi0 = phi;  
    
    if 1 == option.energyflag
        calculate_energy(out1,out2,t,phi0,u0,v0);
    end

    if  0 == mod(nt,nsave)
        if 1 == option.printflag
            timeElapsed = toc(tstart);
            fprintf('epsilon=%.3f,t=%.4f/%.3f, dt=%.2e, Nx=%d, Ny=%d, timeElapsed=%f\n',epsilon,t,T,dt,Nx,Ny,timeElapsed);
        end
        
        if 1 == option.saveflag
            ss = [dir_data '/phi_t=' num2str(t) '.txt'];
            fid = fopen(ss, 'wt');
            fprintf(fid, '%f\n', phi(:));
            fclose(fid);
        end
        
        nfigure = nfigure +1;
        if 1 == option.plotflag
            if 1 == option.vtkflag
                write_vtk_grid_values(dir_data,x,y,nt,phi0);
            end
            if 1 == option.saveflag
                showsolution_2D(nfigure,xx,yy,phi,t,dir_fig);
            else
%                 subplot(1,2,1);
                showsolution_2D(nfigure,xx,yy,phi,t);
%                 subplot(1,2,2);
%                 showsolution_2D(nfigure,xx,yy,pde.exact(xx,yy,t),t);
            end
        end
    end
    
end

if 1 == option.savefinal
    name=['phi_e',num2str(pde.epsilon),'M',num2str(pde.M),'S1=',num2str(pde.S1),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];    
    filename=[name '.mat'];
    save(filename,'epsilon','xx','yy','hx','hy','Nx','Ny','dt','T','phi','domain');
end

if 1 == option.energyflag
    fclose(out1);
    fclose(out2);
end
end

function result = fun_u_init(phi)
global epsilon C0
E1 = fun_inner(1,(lap_diff(phi)-1./epsilon.^2*f(phi)).^2);
if E1 + C0 <0
    disp("Root < 0");
    return;
end
result  = sqrt(E1 + C0);
end

function result = fun_v_init(phi)
global beta
result  = B(phi) - beta;
end

function result = fun_H(phi)
global epsilon C0 Lx Ly
E1 = fun_inner(1,(lap_diff(phi)-1./epsilon.^2*f(phi)).^2);
if E1 + C0 < 0
    disp("Root < 0");
    return;
end
delta_V = -lap_diff(omega(phi)) + 1/epsilon.^2*f_der(phi).*omega(phi);
r1 = delta_V./sqrt(E1 + C0);
r2 = fun_inner(1,r1);
result  = r1 - r2/(Lx*Ly);
end

function result = fun_G(phi)
global epsilon Lx Ly
r1 = -lap_diff(phi) + (1/epsilon.^2).*f(phi);
r2 = fun_inner(1,r1);
result  = r1 - r2/(Lx*Ly);
end

function result = fun_u(phi,phi0,u0,H)
Hphi0 = fun_inner(H,phi0);
Hphi1 = fun_inner(H,phi);
g1 = u0 - Hphi0;
result  = Hphi1 + g1;
end

function result = fun_v(phi,phi0,v0,G)
global epsilon
Gphi0 = fun_inner(G,phi0);
Gphi1 = fun_inner(G,phi);
g2 = v0 - epsilon.*Gphi0;
result  = epsilon.*Gphi1 + g2;
end

function result = get_rhs(phi0,u0,v0,H,G)
global dt M epsilon M2 S1 S2 S3
Hphi0 = fun_inner(H,phi0);
Gphi0 = fun_inner(G,phi0);
g1 = u0 - Hphi0;
g2 = v0 - epsilon.*Gphi0;

result = phi0/epsilon/dt/M - H.*g1 - M2.*G.*g2 ...
         + S1/epsilon.^4.*phi0 - S2/epsilon.^2*lap_diff(phi0) ...
         + S3*lap_diff(lap_diff(phi0));
end

function result = inv_A(phi)
global dt k2 M epsilon S1 S2 S3
    L1 = S1/epsilon.^4 - S2/epsilon.^2*k2 + S3*k2.^2;
    phihat = fft2(phi);
    r      = phihat./(1/epsilon/dt/M + L1);
    result = real(ifft2(r));
end

function r = fun_inner(f,g)
global hx hy
    r1 = fft2(f.*g);
    r = r1(1,1)*hx*hy;
end

function [] = calculate_energy(out1,out2,t,phi,u,v)
global C0 epsilon M2 beta

energy_part1 = fun_inner(1,epsilon./2.*(lap_diff(phi) - f(phi)).^2) ;

energy_part2 = 1./2.*M2.*(B(phi) - beta).^2;

energy_original = energy_part1 + energy_part2;

energy_modified = epsilon./2.*(u.^2 - C0) + 1./2.*M2.*v.^2;

mass    = fun_inner(1,phi);
surface = B(phi);

fprintf(out1,'%14.6e  %.10f %.10f \n',t,mass,surface);
fprintf(out2,'%14.6e  %f  %f\n',t,energy_original,energy_modified);
end

function lap=lap_diff(phi)
global k2
    lap=real(ifft2((k2.*fft2(phi))));
end

function lap=diff_x(phi)
global kx
    lap=real(ifft2((kx.*fft2(phi))));
end

function lap=diff_y(phi)
global ky
    lap=real(ifft2((ky.*fft2(phi))));
end

function result = grad_square(phi)
    result = diff_x(phi).^2 + diff_y(phi).^2;
end

function result = omega(phi)
global epsilon
    result = -lap_diff(phi) + 1/(epsilon.^2).*f(phi);
end

function result = f_der(phi)
    result = 3.*phi.^2-1;
end

function result = f(phi)
    result = phi.^3 - phi;
end

function result = F(phi)
    result = 1/4*(phi.^2-1).^2;
end

function result = B(phi)
global epsilon
    result = fun_inner(1,epsilon./2.*grad_square(phi)+1/epsilon.*F(phi));
end

