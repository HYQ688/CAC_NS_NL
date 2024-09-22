function [phi, u, v, p, Q, eta_1, eta_2] = CAC_Vesicle_with_NS_2D_newLM_p_SAV_bdf2(pde,domain,Nx,Ny,time,option)
global dt kx ky kxx kyy k2 k4 hx hy Lx Ly ...
       epsilon gamma Re lambda C0 S1 S2 S3

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
    figname_mass = [pde.name, '_dt', num2str(time.dt), '_S1', num2str(pde.S1), '_mass_data.txt'];
    figname_energy = [pde.name, '_dt', num2str(time.dt), '_S1', num2str(pde.S1),'_energy.txt'];      
    out1 = fopen(figname_mass,'a');
    out2 = fopen(figname_energy,'a');
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
gamma   = pde.gamma;
lambda  = pde.lambda;
Re      = pde.Re;
C0      = pde.C0;
S1      = pde.S1;
S2      = pde.S2;
S3      = pde.S3;

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
phi0 = pde.initphi(xx,yy);
phi_init = phi0;
u0 = pde.initu(xx,yy);
v0 = pde.initv(xx,yy);
p0 = pde.initp(xx,yy);
nfigure =1;

nplot = round((T-t)/dt);
nsave = round(tsave/dt);

tstart = tic;

mu0 = lambda.*( epsilon.*lap_diff(lap_diff(phi0)) + w(phi0));
Q0  = 1;

%% initialization phi 1
time1 = time; time1.T = dt;
[phi1,u1,v1,p1,Q1,~,~] = CAC_Vesicle_with_NS_2D_newLM_p_SAV_1st(pde,domain,Nx,Ny,time1,option);
mu1 = lambda.*( epsilon.*lap_diff(lap_diff(phi1)) + w(phi1));
t = t+dt;

for nt = 2:nplot
    t = t+dt;
     
    phi_star = 2*phi1-phi0;
    mu_star = 2*mu1-mu0;
    u_star =  2*u1-u0;
    v_star =  2*v1-v0;
    
    % step 1
    if isfield(pde,'rhs') && isfield(pde,'exact')
%         ephi   = pde.exact(xx,yy,t);
%         ephi_t = pde.exact_t(xx,yy,t);
%         emu   = epsilon*lap_diff(lap_diff(ephi)) + fun_q(ephi);
%         tmp  = ephi_t./M - lap_diff(emu);
        rhs = pde.rhs(xx,yy,t);
%         rhs = rhs./M;
    else
        rhs = 0;
    end

    g11 = -(delta_W(phi_star)-1./(Lx*Ly).*fun_inner(delta_W(phi_star),1));
    g12 = -(u_star.*diff_x(phi_star) + v_star.*diff_y(phi_star))./gamma./lambda;
    g21 = (delta_A(phi_star)-1./(Lx*Ly).*fun_inner(delta_A(phi_star),1));
    g22 = g12;
    g31 = (4*phi1-phi0)/2/dt/gamma/lambda + rhs...
           + S1./epsilon.^3.*(phi_star - 1./(Lx.*Ly).*fun_inner(phi_star,1)) ...
           - S2./epsilon.*(lap_diff(phi_star) - 1./(Lx.*Ly).*fun_inner(lap_diff(phi_star),1)) ...
           + S3.*epsilon.*(lap_diff(lap_diff(phi_star)) - 1./(Lx.*Ly).*fun_inner(lap_diff(lap_diff(phi_star)),1));
    g32 = g12;
    
    phi_11 = inv_A(g11);
    phi_12 = inv_A(g12);
    phi_21 = inv_A(g21);
    phi_22 = inv_A(g22);
    phi_31 = inv_A(g31);
    phi_32 = inv_A(g32);

    mu_11 = get_mu1(phi_11) + lambda.*delta_W(phi_star);
    mu_12 = get_mu1(phi_12);
    mu_21 = get_mu1(phi_21) - lambda.*delta_A(phi_star);
    mu_22 = get_mu1(phi_22);
    mu_31 = get_mu2(phi_31,phi_star);
    mu_32 = get_mu1(phi_32);

    %step 2
    G11_u = 0; G11_v = 0;
    G12_u = -(u_star.*diff_x(u_star) + v_star.*diff_y(u_star)) + mu_star.*diff_x(phi_star);
    G12_v = -(u_star.*diff_x(v_star) + v_star.*diff_y(v_star)) + mu_star.*diff_y(phi_star);
    G21_u = 0; G21_v = 0;
    G22_u = G12_u;       G22_v = G12_v;
    G31_u = (4*u1-u0)/2./dt-diff_x(p1); G31_v = (4*v1-v0)/2./dt-diff_y(p1);
    G32_u = G12_u;       G32_v = G12_v;
    utilde_11 = inv_B(G11_u); vtilde_11 = inv_B(G11_v);
    utilde_12 = inv_B(G12_u); vtilde_12 = inv_B(G12_v);
    utilde_21 = inv_B(G21_u); vtilde_21 = inv_B(G21_v);
    utilde_22 = inv_B(G22_u); vtilde_22 = inv_B(G22_v);
    utilde_31 = inv_B(G31_u); vtilde_31 = inv_B(G31_v);
    utilde_32 = inv_B(G32_u); vtilde_32 = inv_B(G32_v);

    %step 3
    nu_11 = get_nu(utilde_11,vtilde_11,mu_11,phi_star,mu_star,u_star,v_star);
    nu_12 = get_nu(utilde_12,vtilde_12,mu_12,phi_star,mu_star,u_star,v_star);
    Q_1 = nu_11/(3-nu_12);
    nu_21 = get_nu(utilde_21,vtilde_21,mu_21,phi_star,mu_star,u_star,v_star);
    nu_22 = get_nu(utilde_22,vtilde_22,mu_22,phi_star,mu_star,u_star,v_star);
    Q_2 = nu_21/(3-nu_22);
    nu_31 = get_nu(utilde_31,vtilde_31,mu_31,phi_star,mu_star,u_star,v_star) + (4*Q1-Q0);
    nu_32 = get_nu(utilde_32,vtilde_32,mu_32,phi_star,mu_star,u_star,v_star);
    Q_3 = nu_31/(3-nu_32);

    %step 4
    phi_1 = phi_11 + Q_1*phi_12; phi_2 = phi_21 + Q_2*phi_22; phi_3 = phi_31 + Q_3*phi_32;
     mu_1 = mu_11 + Q_1*mu_12;    mu_2 = mu_21 + Q_2*mu_22;    mu_3 = mu_31 + Q_3*mu_32;
    utilde_1 = utilde_11 + Q_1*utilde_12; utilde_2 = utilde_21 + Q_2*utilde_22; utilde_3 = utilde_31 + Q_3*utilde_32;
    vtilde_1 = vtilde_11 + Q_1*vtilde_12; vtilde_2 = vtilde_21 + Q_2*vtilde_22; vtilde_3 = vtilde_31 + Q_3*vtilde_32;
    
    %step 5
    eta_2_initial = fun_inner(delta_A(phi_star),4*phi1-phi0-3*(phi_1+phi_3)) ./ fun_inner(delta_A(phi_star),3*phi_2);
    % eta_2_initial = fun_inner(delta_A(phi1),phi1-phi_1-phi_3) ./ fun_inner(delta_A(phi1),phi_2);
    eta_2 = fsolve(@(eta_2)non_fun1(eta_2,phi_1,phi_2,phi_3,phi_init),eta_2_initial,...
                   opts);
    
    %step 6
    eta_1_initial = 1;
    eta_1 = fsolve(@(eta_1)non_fun2(eta_1,eta_2,phi_1,phi_2,phi_3,phi1,phi0),eta_1_initial,...
                   opts);
   
    %step 7
    phi = eta_1.*phi_1+eta_2.*phi_2+phi_3;
    mu  = eta_1.*mu_1+eta_2.*mu_2+mu_3;
    Q  = eta_1.*Q_1+eta_2.*Q_2+Q_3;
    utilde = eta_1.*utilde_1+eta_2.*utilde_2+utilde_3;
    vtilde = eta_1.*vtilde_1+eta_2.*vtilde_2+vtilde_3;

    %step 8
    p1_hat = fft2(p1);
    ut_hat = fft2(utilde);
    vt_hat = fft2(vtilde);
    p_hat = p1_hat + (kx.*ut_hat+ky.*vt_hat)*3/(2*dt)./k2;
    p_hat(1,1) = 0;
    p = real(ifft2(p_hat));

    u = utilde -2*dt/3*diff_x(p-p1);
    v = vtilde -2*dt/3*diff_y(p-p1);

    %% update phi0
    phi0 = phi1; phi1 = phi;
    mu0  = mu1;  mu1  = mu;
    u0   = u1;   u1   = u;
    v0   = v1;   v1   = v;
    p0   = p1;   p1   = p;
    Q0   = Q1;   Q1   = Q;

    if 1 == option.energyflag
        calculate_energy1(out1,out2,t,phi1,phi0,u1,u0,v1,v0,p1,Q1,Q0);
    end

    if  0 == mod(nt,nsave)
        if 1 == option.printflag
            timeElapsed = toc(tstart);
            fprintf('Q=%.4e,eta_1=%.4e,eta_2=%.4e,epsilon=%.3f,t=%.5f/%.4f, dt=%.2e, Nx=%d, Ny=%d, timeElapsed=%f\n',Q,eta_1,eta_2,epsilon,t,T,dt,Nx,Ny,timeElapsed);
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
                write_vtk_grid_values(dir_data,x,y,nt,phi0);
            end
            if 1 == option.saveflag
                showsolution_2D(nfigure,xx,yy,phi,t,dir_fig);
            else
                showsolution_2D(nfigure,xx,yy,phi,t);
            end
        end
    end
    
end

if 1 == option.savefinal
    name=[option.scheme,'_phi_e',num2str(pde.epsilon),'gamma',num2str(pde.gamma),'S1',num2str(pde.S1),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];    
    filename=[name '.mat'];
    save(filename,'epsilon','xx','yy','hx','hy','Nx','Ny','dt','T','phi','u','v','p','domain');
end

if 1 == option.energyflag
    fclose(out1);
    fclose(out2);
end
end

function result = fun_U_init(phi)
global C0
if fun_inner(1,W(phi)) + C0 <0
    disp("Root < 0");
    return;
end
result  = sqrt(fun_inner(1,W(phi)) + C0);
end


function result = inv_A(phi)
global dt k2 gamma lambda epsilon Lx Ly S1 S2 S3
    L1 = epsilon.*k2.^2;
    L2 = S1/epsilon.^3 - S2/epsilon.*k2 + S3.*epsilon.*k2.^2;
    phihat = fft2(phi);
    r      = phihat./(3/2/dt/gamma/lambda + L1 + L2);
    r(1,1) = phihat(1,1)./(3/2/dt/gamma/lambda + L1(1,1) + L2(1,1) - 1/(Lx*Ly).*L1(1,1).*Lx.*Ly - 1/(Lx*Ly).*L2(1,1).*Lx.*Ly);
    result = real(ifft2(r));
end

function result = inv_B(u)
global dt Re k2
    uhat = fft2(u);
    r = uhat./(3/2/dt - 1/Re.*k2);
    result = real(ifft2(r));
end

function result = W(phi)
global epsilon 
    result = 6./epsilon.^2.*phi.^2.*grad_square(phi) ...
             - 2/epsilon.^2.*grad_square(phi) ...
             + 1/epsilon.^4.*(f(phi)).^2;
    result = result*epsilon/2; 
    result = fun_inner(1,result);
end

function result = w(phi)
global epsilon
    div_term = diff_x(phi.^2.*diff_x(phi)) + diff_y(phi.^2.*diff_y(phi));
    result = 6./epsilon.*(phi.*grad_square(phi) - div_term) ...
             + 2/epsilon.*lap_diff(phi) ...
             + 1/epsilon.^3.*f(phi).*f_der(phi);
end

function result = delta_W(phi)
global epsilon
    div_term = diff_x(phi.^2.*diff_x(phi)) + diff_y(phi.^2.*diff_y(phi));
    result = 6./epsilon.*(phi.*grad_square(phi) - div_term) ...
             + 2/epsilon.*lap_diff(phi) ...
             + 1/epsilon.^3.*f(phi).*f_der(phi);
end

function result = get_mu1(phi)
    global epsilon lambda S1 S2 S3
    result = lambda.*( epsilon.*lap_diff(lap_diff(phi)) + S1./epsilon.^3.*phi - S2./epsilon.*lap_diff(phi) + epsilon.*S3.*lap_diff(lap_diff(phi))  );
end

function result = get_mu2(phi,phi_star)
    global epsilon lambda S1 S2 S3
    result = lambda.*( epsilon.*lap_diff(lap_diff(phi)) + S1./epsilon.^3.*(phi - phi_star) - S2./epsilon.*lap_diff(phi - phi_star) + epsilon.*S3.*lap_diff(lap_diff(phi - phi_star))  );
end

function result = get_nu(u,v,mu,phi_star,mu_star,u_star,v_star)
global dt hx hy
    result = fft2((u_star.*diff_x(phi_star)+v_star.*diff_y(phi_star)).*mu ...
                 -mu_star.*((diff_x(phi_star).*u + diff_y(phi_star).*v)) ...
                 +(u_star.*diff_x(u_star) + v_star.*diff_y(u_star)).*u ...
                 +(u_star.*diff_x(v_star) + v_star.*diff_y(v_star)).*v) ;
    result = 2*dt.*result(1,1)*hx*hy;
end

function r = fun_inner(f,g)
global hx hy
    r1 = fft2(f.*g);
    r = r1(1,1)*hx*hy;
end

function [] = calculate_energy1(out1,out2,t,phi1,phi0,u1,u0,v1,v0,p1,Q1,Q0)
global epsilon lambda S1 S2 S3 dt

energy_part1 = lambda*fun_inner(1,1./2.*epsilon.*lap_diff(phi1).^2);

energy_part2 = 1./2*fun_inner(1,u1.^2 + v1.^2);

energy_part3 = lambda*W(phi1);

energy_original = energy_part1 + energy_part2 + energy_part3;

energy_phi = fun_inner(1,lap_diff(phi1)^2) + fun_inner(1,(2*lap_diff(phi1)-lap_diff(phi0))^2);
energy_phi = lambda*epsilon*energy_phi./4;
% energy_phi = fun_inner(1,epsilon/4.*(lap_diff(phi1).^2 + lap_diff(phi_star).^2));


energy_u = fun_inner(1,u1.^2 + v1.^2) + fun_inner(1,((2*u1-u0)^2+(2*v1-v0)^2));
energy_u = energy_u./4;

energy_p = fun_inner(1,dt.^2/3*grad_square(p1));

energy_W = lambda*(3*W(phi1) - W(phi0))./2;

energy_Q = 1./2*(0.5*Q1.^2 + 0.5*(2*Q1-Q0).^2) - 0.5;

energy_S = S1./2./epsilon.^3.*(phi1 - phi0).^2 ...
         + S2./2./epsilon.*grad_square(phi1 - phi0) ...
         + S3.*epsilon./2.*lap_diff(phi1 - phi0).^2;
energy_S = lambda*fun_inner(1,energy_S);

energy_original_discrete = energy_part1 + energy_u + energy_p + energy_W + energy_Q + energy_S;

mass    = fun_inner(1,(phi1 + 1)./2);

surface = A(phi1);

fprintf(out1,'%14.6e  %.8f %.8f \n',t,mass,surface);
fprintf(out2,'%14.6e  %f  %f \n',t,energy_original,energy_original_discrete);
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

function result = f_der(phi)
    result = 3.*phi.^2-1;
end

function result = f(phi)
    result = phi.^3 - phi;
end

function result = F(phi)
    result = 1/4*(phi.^2-1).^2;
end

function result = A(phi)
global hx hy epsilon
    r = fft2(epsilon./2.*grad_square(phi)+1/epsilon.*F(phi));
    result = r(1,1)*hx*hy;
end

function r = b(phi)
global epsilon
    r = epsilon./2.*grad_square(phi)+1/epsilon.*F(phi);
end

function result = delta_A(phi)
global epsilon
    result = -epsilon*lap_diff(phi) + 1/epsilon.*f(phi);
end

function r = non_fun1(eta_2,phi_1,phi_2,phi_3,phi_init)
global hx hy
    left = fft2(b(phi_1 + eta_2.*phi_2 + phi_3));
    right = fft2(b(phi_init));
    r = left(1,1)*hx*hy - right(1,1)*hx*hy ;
end

function r = non_fun2(eta_1,eta_2,phi_1,phi_2,phi_3,phi1,phi0)
global hx hy
    phi_star = 2*phi1 - phi0 ;
    left = fft2(3.*W(eta_1.*phi_1+eta_2.*phi_2+phi_3)-4.*W(phi1)+W(phi0));
%     left = 3.*W(phi_1+eta.*phi_2+lambda.*phi_3) - 4.*W(phi1) + W(phi0); 
    right1 = eta_1.*fft2(delta_W(phi_star).*(3.*(eta_1.*phi_1+eta_2.*phi_2+phi_3)-4.*phi1+phi0));
    right2 = eta_2.*fft2(delta_A(phi_star).*(3.*(eta_1.*phi_1+eta_2.*phi_2+phi_3)-4.*phi1+phi0));
    r = left(1,1)*hx*hy - right1(1,1)*hx*hy  + right2(1,1)*hx*hy ;
end



