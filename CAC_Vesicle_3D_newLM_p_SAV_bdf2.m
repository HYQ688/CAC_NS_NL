function [phi, eta, lambda] = CAC_Vesicle_3D_newLM_p_SAV_bdf2(pde,domain,Nx,Ny,Nz,time,option)

global dt kx ky kz kxx kyy kzz k2 k4 hx hy hz Lx Ly Lz ...
       epsilon M C0 S1 S2 S3

if ~exist('option','var'), option = []; end
if ~isfield(option,'tol')
    option.tol = 10^-14;   % default tol
end
if ~isfield(option,'tolit')
    option.tolit = 10^-11;   % default tolit
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

%% options for fsolve
opts = optimoptions('fsolve','Display','off',...
                       'StepTolerance', option.tol, ...
                       'FunctionTolerance',option.tolit,...
                       'MaxIterations',option.maxit);

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
S1      = pde.S1;
S2      = pde.S2;
S3      = pde.S3;

Lx = domain.xb - domain.xa;
Ly = domain.yb - domain.ya;
Lz = domain.zb - domain.za;

hx = Lx/Nx;
hy = Ly/Ny;
hz = Lz/Nz;
% x  = domain.left   + hx*(1:Nx);
% y  = domain.bottom + hy*(1:Ny);
x = linspace( domain.xa, domain.xb-hx, Nx); 
y = linspace( domain.ya, domain.yb-hy, Ny);
z = linspace( domain.za, domain.zb-hz, Nz);

% [k_x,k_y,kx,ky,kxx,kyy,k2,k4] = prepare_fft3_v2(Lx,Ly,Lz,Nx,Ny,Nz);
k_x = 1i*[0:Nx/2 -Nx/2+1:-1]*(2*pi/Lx);
k_y = 1i*[0:Ny/2 -Ny/2+1:-1]*(2*pi/Ly);
k_z = 1i*[0:Nz/2 -Nz/2+1:-1]*(2*pi/Lz);
[kx, ky, kz] = ndgrid(k_x,k_y,k_z);

k2x = k_x.^2;
k2y = k_y.^2;
k2z = k_z.^2;
[kxx, kyy, kzz] = ndgrid(k2x,k2y,k2z);
k2 = kxx + kyy + kzz;
k4 = k2.^2;

[xx,yy,zz] = ndgrid(x,y,z);
phi0 = pde.init(xx,yy,zz);      
phi_init = phi0;
nfigure =1;

nplot = round((T-t)/dt);
nsave = round(tsave/dt);

tstart = tic;

%% initialization phi 1
time1 = time; time1.T = dt;
[phi1,~,lambda0] = CAC_Vesicle_3D_newLM_p_SAV_1st(pde,domain,Nx,Ny,Nz,time1,option);

t = t+dt;

for nt = 2:nplot
    t = t+dt;
     
    phi_star = 2*phi1-phi0;    
    
    % step 1
    rhs1 = (4*phi1-phi0)/2/dt/M...
           + S1./epsilon.^3.*(phi_star - 1./(Lx.*Ly.*Lz).*fun_inner(phi_star,1)) ...
           - S2./epsilon.*(lap_diff(phi_star) - 1./(Lx.*Ly.*Lz).*fun_inner(lap_diff(phi_star),1)) ...
           + S3.*epsilon.*(lap_diff(lap_diff(phi_star)) - 1./(Lx.*Ly.*Lz).*fun_inner(lap_diff(lap_diff(phi_star)),1));
    
    if isfield(pde,'rhs') && isfield(pde,'exact')
%         ephi   = pde.exact(xx,yy,t);
%         ephi_t = pde.exact_t(xx,yy,t);
%         emu   = epsilon*lap_diff(lap_diff(ephi)) + fun_q(ephi);
%         tmp  = ephi_t./M - lap_diff(emu);
        rhs = pde.rhs(xx,yy,zz,t);
%         rhs = rhs./M;
    else
        rhs = 0;
    end

    g1 = rhs1 + rhs;
    g2 = -(delta_W(phi_star)-1./(Lx*Ly*Lz).*fun_inner(delta_W(phi_star),1));
    g3 = (delta_B(phi_star)-1./(Lx*Ly*Lz).*fun_inner(delta_B(phi_star),1));

    phi_1 = inv_A(g1);
    phi_2 = inv_A(g2);
    phi_3 = inv_A(g3);

    %step 2
    lambda_initial = fun_inner(delta_B(phi_star),4*phi1-phi0-3*(phi_1+phi_2)) ./ fun_inner(delta_B(phi_star),3*phi_3);
    lambda = fsolve(@(lambda)non_fun1(lambda,phi_1,phi_2,phi_3,phi_init),lambda_initial,...
                   opts);

    %step 3
    eta_initial = 1;
    eta = fsolve(@(eta)non_fun2(eta,lambda,phi_1,phi_2,phi_3,phi1,phi0),eta_initial,...
                   opts);
%     eta = 1;

    %step 4
    phi = phi_1+eta.*phi_2+lambda.*phi_3;

    %% update phi0
    phi0 = phi1; 
    phi1 = phi;
    
    if 1 == option.energyflag
        calculate_energy1(out1,out2,t,phi1,phi0);
    end

    if  0 == mod(nt,nsave)
        if 1 == option.printflag
            timeElapsed = toc(tstart);
            fprintf('eta=%.4e,lambda=%.4e,epsilon=%.3f,t=%.5f/%.4f, dt=%.2e, Nx=%d, Ny=%d, Nz=%d,timeElapsed=%f\n',eta,lambda,epsilon,t,T,dt,Nx,Ny,Nz,timeElapsed);
        end
        
        if 1 == option.saveflag
            if ~exist(dir_data,'dir')
                mkdir(dir_data);
            end
            ss1 = [dir_data '/phi_t=' num2str(t) '.txt'];
            writematrix(phi0,ss1,'Delimiter',' ');
        end
        
        nfigure = nfigure +1;
        if 1 == option.plotflag
            if 1 == option.vtkflag
                write_vtk_grid_values(dir_data,x,y,z,nt,phi0);
            end
            if 1 == option.saveflag
                showsolution_3D(nfigure,xx,yy,zz,phi,t,dir_fig);
            else
                showsolution_3D(nfigure,xx,yy,zz,phi,t);
            end
        end
    end
    
end

if 1 == option.savefinal
    name=[option.scheme,'_phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'dt=',num2str(dt)];    
    filename=[name '.mat'];
    save(filename,'epsilon','xx','yy','zz','hx','hy','hz','Nx','Ny','Nz','dt','T','phi','domain');
end

if 1 == option.energyflag
    fclose(out1);
    fclose(out2);
end
end

function result = fun_U_init(phi)
global C0
if fun_inner(1,fun_W(phi)) + C0 <0
    disp("Root < 0");
    return;
end
result  = sqrt(fun_inner(1,fun_W(phi)) + C0);
end


function result = inv_A(phi)
global dt k2 M epsilon Lx Ly Lz S1 S2 S3
    L1 = epsilon.*k2.^2;
    L2 = S1/epsilon.^3 - S2/epsilon*k2 + S3.*epsilon.*k2.^2;
    phihat = fftn(phi);
    r      = phihat./(3/2/dt/M + L1 + L2);
    r(1,1,1) = phihat(1,1,1)./(3/2/dt/M);
    result = real(ifftn(r));
end

function result = fun_W(phi)
global epsilon 
    result = 6./epsilon.^2.*phi.^2.*grad_square(phi) ...
             - 2/epsilon.^2.*grad_square(phi) ...
             + 1/epsilon.^4.*(f(phi)).^2;
    result = result*epsilon/2;    
end

function result = delta_W(phi)
global epsilon
    div_term = diff_x(phi.^2.*diff_x(phi)) + diff_y(phi.^2.*diff_y(phi)) + diff_z(phi.^2.*diff_z(phi));
    result = 6./epsilon.*(phi.*grad_square(phi) - div_term) ...
             + 2/epsilon.*lap_diff(phi) ...
             + 1/epsilon.^3.*f(phi).*f_der(phi);
end

function r = fun_inner(f,g)
global hx hy hz
    r1 = fftn(f.*g);
    r = r1(1,1,1)*hx*hy*hz;
end

function [] = calculate_energy1(out1,out2,t,phi1,phi0)
global epsilon S1 S2 S3

phi_star = 2*phi1 - phi0;
 
energy_linear = fun_inner(1,1./2.*epsilon.*lap_diff(phi1).^2);

energy_nonlinear = fun_inner(1,fun_W(phi1));

energy_original = energy_linear + energy_nonlinear;

energy_S = S1./2./epsilon.^3.*(phi1 - phi0).^2 ...
         + S2./2./epsilon.*grad_square(phi1 - phi0) ...
         + S3.*epsilon./2.*lap_diff(phi1 - phi0).^2;
energy_S = fun_inner(1,energy_S);

energy_part1 = fun_inner(1,epsilon/4.*(lap_diff(phi1).^2 + lap_diff(phi_star).^2));
energy_part2 = (3*fun_inner(1,fun_W(phi1)) - fun_inner(1,fun_W(phi0)))./2;

energy_original_discrete = energy_part1 + energy_part2 + energy_S;

mass    = fun_inner(1,(phi1 + 1)./2);

surface = B(phi1);

fprintf(out1,'%14.6e  %.8f %.8f \n',t,mass,surface);
fprintf(out2,'%14.6e  %f  %f \n',t,energy_original,energy_original_discrete);
end

function lap=lap_diff(phi)
global k2
    lap=real(ifftn((k2.*fftn(phi))));
end

function lap=diff_x(phi)
global kx
    lap=real(ifftn((kx.*fftn(phi))));
end

function lap=diff_y(phi)
global ky
    lap=real(ifftn((ky.*fftn(phi))));
end

function lap=diff_z(phi)
global kz
    lap=real(ifftn((kz.*fftn(phi))));
end

function result = grad_square(phi)
    result = diff_x(phi).^2 + diff_y(phi).^2 + diff_z(phi).^2;
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

function r = b(phi)
global epsilon
    r = epsilon./2.*grad_square(phi)+1/epsilon.*F(phi);
end

function result = delta_B(phi)
global epsilon
    result = -epsilon*lap_diff(phi) + 1/epsilon.*f(phi);
end

function r = W(phi)
global hx hy hz
    r = fftn(fun_W(phi));
    r = r(1,1,1)*hx*hy*hz;
end

function r = non_fun1(lambda,phi_1,phi_2,phi_3,phi0)
global hx hy hz
    left = fftn(b(phi_1 + phi_2 + lambda.*phi_3));
    right = fftn(b(phi0));
    r = left(1,1)*hx*hy*hz - right(1,1)*hx*hy*hz ;
end

function r = non_fun2(eta,lambda,phi_1,phi_2,phi_3,phi1,phi0)
global hx hy hz
    phi_star = 2*phi1 - phi0 ;
    left = fftn(3.*fun_W(phi_1+eta.*phi_2+lambda.*phi_3)-4.*fun_W(phi1)+fun_W(phi0));
%     left = 3.*W(phi_1+eta.*phi_2+lambda.*phi_3) - 4.*W(phi1) + W(phi0); 
    right1 = eta.*fftn(delta_W(phi_star).*(3.*(phi_1+eta.*phi_2+lambda.*phi_3)-4.*phi1+phi0));
    right2 = lambda.*fftn(delta_B(phi_star).*(3.*(phi_1+eta.*phi_2+lambda.*phi_3)-4.*phi1+phi0));
    r = left(1,1)*hx*hy*hz - right1(1,1)*hx*hy*hz  + right2(1,1)*hx*hy*hz ;
end



