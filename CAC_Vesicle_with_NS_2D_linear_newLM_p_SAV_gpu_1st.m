function [phi, U, u, v, p, Q, eta_2] = CAC_Vesicle_with_NS_2D_linear_newLM_p_SAV_gpu_1st(pde,domain,Nx,Ny,time,option)

global dt kx ky kxx kyy k2 k4 hx hy Lx Ly ...
       epsilon gamma Re lambda C0 S1 S2 S3

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
u0 = pde.initu(xx,yy);
v0 = pde.initv(xx,yy);
p0 = pde.initp(xx,yy);
phi0 = gpuArray(phi0);
u0 = gpuArray(u0);
v0 = gpuArray(v0);
p0 = gpuArray(p0);
nfigure =1;

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

mu0 = lambda.*( epsilon.*lap_diff(lap_diff(phi0)) + w(phi0));
U0 = U_init(phi0);
Q0  = 1;

% Initial energy
if 1 == option.energyflag
    calculate_energy1(out1,out2,t,u0,v0,phi0,U0,Q0);
end

for nt = 1:nplot
    t = t+dt;
     
    phi_star = phi0;
    mu_star = mu0;
    u_star = u0;
    v_star = v0;
    
    % step 1
    H = fun_H(phi_star);
    K = U0 - 1/2*fun_inner(H,phi0);

    %step 2 

    g1 = phi0./dt./gamma./lambda ...
            + S1./epsilon.^3.*(phi_star - 1./(Lx.*Ly).*fun_inner(phi_star,1)) ...
            - S2./epsilon.*(lap_diff(phi_star) - 1./(Lx.*Ly).*fun_inner(lap_diff(phi_star),1)) ...
            + S3.*epsilon.*(lap_diff(lap_diff(phi_star)) - 1./(Lx.*Ly).*fun_inner(lap_diff(lap_diff(phi_star)),1));
    g2 = -(u_star.*diff_x(phi_star) + v_star.*diff_y(phi_star))./gamma./lambda;
    g3 = delta_B(phi_star) - 1./(Lx.*Ly).*fun_inner(1,delta_B(phi_star));
    g4 = g2;

    psi1 = inv_A(g1-H.*K); psi2 = inv_A(g2);
    psi3 = inv_A(g3);      psi4 = inv_A(g4);

    psih = inv_A(1./2.*H);

    b1 = fun_inner(psi1,H); b2 = fun_inner(psi2,H);
    b3 = fun_inner(psi3,H); b4 = fun_inner(psi4,H);

    a = fun_inner(psih,H);

    X = b1./(1+a); Y = b2./(1+a);
    Z = b3./(1+a); W = b4./(1+a);

    %step 3
    phi_11 = psi1 - psih.*X; phi_12 = psi2 - psih.*Y;
    phi_21 = psi3 - psih.*Z; phi_22 = psi4 - psih.*W;

    U_11 = K + 1./2.*X; U_12 = 1./2.*Y;
    U_21 = 1./2.*Z;     U_22 = 1./2.*W;

    mu_11 = get_mu1(phi_11,phi_star,U_11,H);
    mu_12 = get_mu2(phi_12,U_12,H);
    mu_21 = get_mu3(phi_21,phi_star,U_21,H);
    mu_22 = get_mu2(phi_22,U_22,H);

    %step 4    
    if isfield(pde,'rhs_u') && isfield(pde,'exact_u')
        rhs1 = pde.rhs_u(xx,yy,t);
        rhs2 = pde.rhs_v(xx,yy,t);
    else
        rhs1 = 0;
        rhs2 = 0;
    end

    G1_u = rhs1 + u0./dt - diff_x(p0);  G1_v = rhs2 + v0./dt - diff_y(p0);
    G2_u = -(u_star.*diff_x(u_star) + v_star.*diff_y(u_star)) + mu_star.*diff_x(phi_star);
    G2_v = -(u_star.*diff_x(v_star) + v_star.*diff_y(v_star)) + mu_star.*diff_y(phi_star);
    G3_u = 0; G3_v = 0;
    G4_u = G2_u; G4_v = G2_v;

    utilde_11 = inv_B(G1_u); vtilde_11 = inv_B(G1_v);
    utilde_12 = inv_B(G2_u); vtilde_12 = inv_B(G2_v);
    utilde_21 = inv_B(G3_u); vtilde_21 = inv_B(G3_v);
    utilde_22 = inv_B(G4_u); vtilde_22 = inv_B(G4_v);

    %step 5
    nu_11 = get_nu(utilde_11,vtilde_11,mu_11,phi_star,mu_star,u_star,v_star) + Q0;
    nu_12 = get_nu(utilde_12,vtilde_12,mu_12,phi_star,mu_star,u_star,v_star);
    Q_1 = nu_11./(1-nu_12);
    
    nu_21 = get_nu(utilde_21,vtilde_21,mu_21,phi_star,mu_star,u_star,v_star);
    nu_22 = get_nu(utilde_22,vtilde_22,mu_22,phi_star,mu_star,u_star,v_star);
    Q_2 = nu_21./(1-nu_22);

    %stpe 6 
    phi_1 = phi_11 + Q_1*phi_12;  phi_2 = phi_21 + Q_2*phi_22;
    mu_1  = mu_11  + Q_1*mu_12;   mu_2  = mu_21  + Q_2*mu_22;
    U_1   = U_11   + Q_1*U_12;    U_2   = U_21   + Q_2*U_22;
    utilde_1 = utilde_11 + Q_1*utilde_12;vtilde_1 = vtilde_11 + Q_1*vtilde_12;
    utilde_2 = utilde_21 + Q_2*utilde_22;vtilde_2 = vtilde_21 + Q_2*vtilde_22;

    %step 7
    eta_2 = fun_inner(delta_B(phi_star),phi0-phi_1)./ fun_inner(delta_B(phi_star),phi_2);
    
    Q  = Q_1 + eta_2.*Q_2;     phi = phi_1 + eta_2.*phi_2;
    mu = mu_1 + eta_2.*mu_2;   U   = U_1 + eta_2.*U_2;
    utilde  = utilde_1 + eta_2.*utilde_2;  vtilde  = vtilde_1 + eta_2.*vtilde_2;

    %step 8
    p0_hat = fft2(p0);
    ut_hat = fft2(utilde);
    vt_hat = fft2(vtilde);
    p_hat = p0_hat + (kx.*ut_hat+ky.*vt_hat)./k2/dt;
    p_hat(1,1) = 0;
    p = real(ifft2(p_hat));

    u = utilde -dt*diff_x(p-p0);
    v = vtilde -dt*diff_y(p-p0);

    %% update 
    phi0 = phi; 
    mu0 = mu;
    U0 = U;
    u0 = u;
    v0 = v;
    p0 = p;
    Q0 = Q;
        
    if 1 == option.energyflag
        calculate_energy1(out1,out2,t,phi,u,v,U,Q);
    end

    if  0 == mod(nt,nsave)
        if 1 == option.printflag
            timeElapsed = toc(tstart);
            fprintf('Q=%.4e,eta_2=%.4e,epsilon=%.3f,t=%.5f/%.4f, dt=%.2e, Nx=%d, Ny=%d, timeElapsed=%f\n',Q,eta_2,epsilon,t,T,dt,Nx,Ny,timeElapsed);
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
   name=[option.scheme,'_phi_e',num2str(pde.epsilon),'gamma=',num2str(pde.gamma),'S1=',num2str(pde.S1),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];    
    filename=[name '.mat'];
    save(filename,'epsilon','xx','yy','hx','hy','Nx','Ny','dt','T','phi','u','v','p','domain');
end

if 1 == option.energyflag
    fclose(out1);
    fclose(out2);
end
end

function result = U_init(phi)
global C0
if W(phi) + C0 <0
    disp("Root < 0");
    return;
end
result  = sqrt(W(phi) + C0);
end

function result = fun_H(phi)
global C0 hx hy Lx Ly
if W(phi) + C0 <0
    disp("Root < 0");
    return;
end
H  = w(phi)./sqrt(W(phi) + C0);
H_bar = fft2(H);
result = H - 1./(Lx*Ly).*H_bar(1,1)*hx*hy;
end

function result = inv_A(phi)
global dt k2 epsilon gamma lambda Lx Ly S1 S2 S3
    L1 = epsilon.*k2.^2;
    L2 = S1/epsilon.^3 - S2/epsilon.*k2 + S3.*epsilon.*k2.^2;
    phihat = fft2(phi);
    r      = phihat./(1/dt/gamma/lambda + L1 + L2);
    r(1,1) = phihat(1,1)./(1/dt/gamma/lambda + L1(1,1) + L2(1,1) - 1/(Lx*Ly).*L1(1,1).*Lx.*Ly - 1/(Lx*Ly).*L2(1,1).*Lx.*Ly);
    result = real(ifft2(r));
end

function result = inv_B(u)
global dt Re k2
    uhat = fft2(u);
    r = uhat./(1/dt - 1/Re.*k2);
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

function result = get_mu1(phi,phi_star,U,H)
    global epsilon lambda S1 S2 S3
    result = lambda.*( epsilon.*lap_diff(lap_diff(phi)) + H.*U + S1./epsilon.^3.*(phi - phi_star) - S2./epsilon.*lap_diff(phi - phi_star) + epsilon.*S3.*lap_diff(lap_diff(phi - phi_star))  );
end

function result = get_mu2(phi,U,H)
    global epsilon lambda S1 S2 S3
    result = lambda.*( epsilon.*lap_diff(lap_diff(phi)) + H.*U + S1./epsilon.^3.*phi - S2./epsilon.*lap_diff(phi) + epsilon.*S3.*lap_diff(lap_diff(phi))  );
end

function result = get_mu3(phi,phi_star,U,H)
    global epsilon lambda S1 S2 S3
    result = lambda.*( epsilon.*lap_diff(lap_diff(phi)) + H.*U - delta_B(phi_star) + S1./epsilon.^3.*phi - S2./epsilon.*lap_diff(phi) + epsilon.*S3.*lap_diff(lap_diff(phi))  );
end

function result = get_nu(u,v,mu,phi_star,mu_star,u_star,v_star)
global dt hx hy
    result = fft2((u_star.*diff_x(phi_star)+v_star.*diff_y(phi_star)).*mu ...
                 -mu_star.*((diff_x(phi_star).*u + diff_y(phi_star).*v)) ...
                 +(u_star.*diff_x(u_star) + v_star.*diff_y(u_star)).*u ...
                 +(u_star.*diff_x(v_star) + v_star.*diff_y(v_star)).*v) ;
    result = dt.*result(1,1)*hx*hy;
end

function r = fun_inner(f,g)
global hx hy
    r1 = fft2(f.*g);
    r = r1(1,1)*hx*hy;
end

function [] = calculate_energy1(out1,out2,t,phi,u,v,U,Q)
global C0 epsilon lambda

energy_part1 = lambda*fun_inner(1,1./2.*epsilon.*lap_diff(phi).^2);

energy_part2 = 1./2*fun_inner(1,u.^2 + v.^2);

energy_part3 = lambda*fun_inner(1,W(phi));

energy_original = energy_part1 + energy_part2 + energy_part3;

energy_phi = lambda*epsilon*fun_inner(lap_diff(phi),lap_diff(phi));


energy_modified = energy_phi + energy_part2 + lambda*(U.^2 - C0) + 1/2*Q;

mass    = fun_inner(1,(phi + 1)./2);

surface = B(phi);

fprintf(out1,'%14.6e  %.8f %.8f \n',t,mass,surface);
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
global hx hy epsilon
    r = fft2(epsilon./2.*grad_square(phi)+1/epsilon.*F(phi));
    result = r(1,1)*hx*hy;
end

function result = delta_B(phi)
global epsilon
    result = -epsilon*lap_diff(phi) + 1/epsilon.*f(phi);
end


