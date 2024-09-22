close all;
clear; clc;

% add path
addpath('../','-begin');

% Space: Domain and N
% domain.left   = -pi;
% domain.right  =  pi;
% domain.bottom = -pi;
% domain.top    =  pi;

domain.left   = 0;
domain.right  = 2*pi;
domain.bottom = 0;
domain.top    = 2*pi;

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;

N  = 64;
Nx = N;
Ny = N;

% scheme1 = 'linear';
scheme1 = 'nonlinear';

% scheme2 = '_1st';   % First-order scheme
% scheme2 = '_2cn';   % Second-order CN scheme
scheme2 = '_bdf2';  % Second-order BDF scheme

scheme = [scheme1 , scheme2];

% Parameters
para.epsilon = 6*pi/128;
para.gamma = 1;
para.Re = 1;
para.lambda = 0.01;
% para.lambda = 1;
para.S1 = 0;
para.S2 = 0;
para.S3 = 0;

PDE = 'data1';

if 1 == strcmp(PDE,'data1')
    T = 2;
    dt_array = 0.1./2.^(1:6)';
    dt_ref = 1e-5;
    para.C0 = 100;
    pde = ex03_1_Vesicles_data(para);
end

if 1 == strcmp(PDE,'data2')
    T = 0.2;
    dt_array = 0.01./2.^(1:6)';
    dt_ref = 1e-5;
    para.C0 = 100;
    pde = ex03_2_Vesicles_data(para);
end

% Time: dt T

t0 = 0;
tsave = 0.2*T;

maxIt = length(dt_array);

%% option
option.scheme = scheme;
option.plotflag  = 0;
option.printflag = 1;
option.vtkflag  = 0;
option.saveflag  = 0;
option.savefinal  = 1;
option.energyflag = 0;
option.tol = 1e-14;
option.tolit = 1e-11;
option.maxit = 2000;

if 1 == strcmp(scheme,'linear_1st')
    solver_fun = @CAC_Vesicle_with_NS_2D_linear_newLM_p_SAV_1st;
elseif 1 == strcmp(scheme,'linear_2cn')
    solver_fun = @CAC_Vesicle_2D_linear_newLM_p_SAV_2cn;
elseif 1 == strcmp(scheme,'linear_bdf2')
    solver_fun = @CAC_Vesicle_with_NS_2D_linear_newLM_p_SAV_bdf2;
elseif 1 == strcmp(scheme,'nonlinear_1st')
    solver_fun = @CAC_Vesicle_with_NS_2D_newLM_p_SAV_1st;
elseif 1 == strcmp(scheme,'nonlinear_2cn')
    solver_fun = @CAC_Vesicle_2D_newLM_p_SAV_2cn;
elseif 1 == strcmp(scheme,'nonlinear_bdf2')
    solver_fun = @CAC_Vesicle_with_NS_2D_newLM_p_SAV_bdf2;
end


%% Run:
% delete *.mat
% if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
%     time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
%     solver_fun(pde,domain,Nx,Ny,time,option);
% end
% for k = 1:maxIt
%     dt = dt_array(k);
%     time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
%     v2 = solver_fun(pde,domain,Nx,Ny,time,option);
% end

%% Compute order of convergence
error_phi=zeros(maxIt,1);
order_phi=zeros(maxIt,1);
error_u=zeros(maxIt,1);
order_u=zeros(maxIt,1);
error_v=zeros(maxIt,1);
order_v=zeros(maxIt,1);
error_p=zeros(maxIt,1);
order_p=zeros(maxIt,1);
if ~isfield(pde,'exactphi') || ~isfield(pde,'rhs1')
    name=[scheme,'_phi_e',num2str(pde.epsilon),'gamma=',num2str(pde.gamma),'S1=',num2str(pde.S1),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt_ref)];
    filename=[name '.mat'];
    load(filename,'phi','u','v','p');
    phi_exact = phi;
    u_exact = u;
    v_exact = v;
    p_exact = p;
    clear phi u v p;
else
    Lx = domain.right - domain.left;
    Ly = domain.top   - domain.bottom;
    hx = Lx/Nx;
    hy = Ly/Ny;
    x  = domain.left   + hx*(0:Nx-1);
    y  = domain.bottom + hy*(0:Ny-1);
    [xx,yy] = ndgrid(x,y);
    phi_exact = pde.exactphi(xx,yy,T);
    u_exact = pde.exactu(xx,yy,T);
    v_exact = pde.exactv(xx,yy,T);
    p_exact = pde.exactp(xx,yy,T);
end

for k = 1:maxIt
    dt = dt_array(k);
    name=[scheme,'_phi_e',num2str(pde.epsilon),'gamma=',num2str(pde.gamma),'S1=',num2str(pde.S1),...
        'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];
    filenamek=[name '.mat'];
    load(filenamek,'phi','u','v','p','hx','hy');
    err_phi    = fft2((phi_exact - phi).^2);
    error_phi(k,1) = sqrt(err_phi(1,1)*hx*hy);   % L2

    err_u = fft2((u_exact - u).^2 + (v_exact - v).^2);
    error_u(k,1) = sqrt(err_u(1,1)*hx*hy);   % L2

    err_p = fft2((p_exact - p).^2);
    error_p(k,1) = sqrt(err_p(1,1)*hx*hy);   % L2

    clear phi u v p;
end
order_phi(2:maxIt) = log(error_phi(1:maxIt-1)./error_phi(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));
order_u(2:maxIt)   = log(error_u(1:maxIt-1)./error_u(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));
order_p(2:maxIt)   = log(error_p(1:maxIt-1)./error_p(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));

%% Display error and order
% fprintf('\nError_phi\n  dt    &   Error_L2   &  Order \n');
% for k = 1:maxIt
%     fprintf('%.5e  &  %.4e  &  %.2f \n',dt_array(k),error_phi(k),order_phi(k));
% end
% 
% fprintf('\nError_u\n  dt    &   Error_L2   &  Order \n');
% for k = 1:maxIt
%     fprintf('%.5e  &  %.4e  &  %.2f \n',dt_array(k),error_u(k),order_u(k));
% end
% fprintf('\n')
% 
% fprintf('Error_p\n  dt    &   Error_L2   &  Order \n');
% for k = 1:maxIt
%     fprintf('%.5e  &  %.4e  &  %.2f \n',dt_array(k),error_p(k),order_p(k));
% end
% fprintf('\n')
fprintf('   dt \t   & phi_Error_L2 &phi_Order &  u_Error_L2  &  u_Order &  p_Error_L2  &  p_Order \n');
for k = 1:maxIt
    fprintf('%.5e  &  %.4e  &  %.4f  &  %.4e  &  %.4f  &  %.4e  &  %.4f \n',dt_array(k),error_phi(k),order_phi(k),error_u(k),order_u(k),error_p(k),order_p(k));
end
X = ['function:',PDE,'method:',scheme];
disp(X)
fprintf('\n')
%% Plot
figure(2)
hh=loglog(dt_array,error_phi,'*-','LineWidth',3,'MarkerSize',15);
hold on;
hh=loglog(dt_array,error_u, 'p-', 'LineWidth',3,'MarkerSize',20);
hh=loglog(dt_array,error_p, 'b-o','LineWidth',3,'MarkerSize',10);
hh=loglog(dt_array,dt_array,'k-.','LineWidth',3,'MarkerSize',15);

tf = strcmp(scheme2,'_1st');
if tf == 1
    hh=loglog(dt_array,dt_array,'k-.','LineWidth',3,'MarkerSize',10);
    hold on
    h=legend('$L^2$-error: $\phi$','$L^2$-error: $\mathbf{u}$','$L^2$-error: $p$','$\mathcal{O}(\delta t)$','Location','SouthEast');
%     h = legend('Linear scheme: $\phi$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
else 
    hh=loglog(dt_array,dt_array.^2,'k-.','LineWidth',3,'MarkerSize',10);
    hold on
    h = legend('nonLinear scheme: $\phi$','$\mathcal{O}(\delta t^2)$','Location','SouthEast');
end
set(h,'interpreter','latex','FontSize',15);
xlabel('Time step $\delta t$','Interpreter','latex');
%  ylabel('Error\_max')
ylabel('$L^2$ error','Interpreter','latex');

ylim([1e-12 5e-1])

grid on;
hold on;

%% Save error and order
name=[scheme,'_phi_e',num2str(para.epsilon),'gamma=',num2str(para.gamma),'S1=',num2str(para.S1),...
      'Nx=',num2str(N),'Ny=',num2str(N)];
% fileID = fopen([name,'.txt'],'w');
% % fprintf(fileID,'%6s\n','%% Results');
% % fprintf(fileID,'%6s\n','% dt	   &   Error_L2	   &  Order');
% % A = [dt_array error];
% % fprintf(fileID,'%.12f   %.4e   \n',A');
% fprintf(fileID,'%.12f     %.4e      %.2f \n',[dt_array,error,order]');
% fclose(fileID);
fileID = fopen([name,'.txt'],'w');
fprintf(fileID,'%.12f   %.4e    %.2f  %.4e   %.2f  %.4e   %.2f\n',...
     [dt_array,error_phi,order_phi,error_u,order_u,error_p,order_p]');
fclose(fileID);

%% results:

%% linear_1st
% ex03_1_Vesicles_data
%    dt 	   & phi_Error_L2 &phi_Order &  u_Error_L2  &  u_Order &  p_Error_L2  &  p_Order 
% 5.00000e-03  &  1.8942e-03  &  0.0000  &  5.1315e-06  &  0.0000  &  5.4491e-04  &  0.0000 
% 2.50000e-03  &  1.0109e-03  &  0.9059  &  2.7832e-06  &  0.8827  &  2.9143e-04  &  0.9028 
% 1.25000e-03  &  5.2155e-04  &  0.9548  &  1.4494e-06  &  0.9413  &  1.5053e-04  &  0.9531 
% 6.25000e-04  &  2.6342e-04  &  0.9855  &  7.3575e-07  &  0.9782  &  7.6074e-05  &  0.9846 
% 3.12500e-04  &  1.3077e-04  &  1.0104  &  3.6621e-07  &  1.0066  &  3.7776e-05  &  1.0099 
% 1.56250e-04  &  6.3515e-05  &  1.0418  &  1.7811e-07  &  1.0399  &  1.8352e-05  &  1.0416 
% 7.81250e-05  &  2.9655e-05  &  1.0988  &  8.3216e-08  &  1.0978  &  8.5690e-06  &  1.0987 
% 3.90625e-05  &  1.2666e-05  &  1.2273  &  3.5554e-08  &  1.2269  &  3.6600e-06  &  1.2273 

% ex03_2_Vesicles_data
% Q=9.9999e-01,eta_2=8.4550e+00,epsilon=0.147,t=0.20000/0.2000, dt=1.56e-04, Nx=64, Ny=64, timeElapsed=15.389783
%    dt 	   & phi_Error_L2 &phi_Order &  u_Error_L2  &  u_Order &  p_Error_L2  &  p_Order 
% 5.00000e-03  &  9.7812e-03  &  0.0000  &  3.5030e-03  &  0.0000  &  1.3579e-02  &  0.0000 
% 2.50000e-03  &  5.2285e-03  &  0.9036  &  1.8684e-03  &  0.9068  &  7.2271e-03  &  0.9099 
% 1.25000e-03  &  2.7002e-03  &  0.9533  &  9.6334e-04  &  0.9557  &  3.7236e-03  &  0.9567 
% 6.25000e-04  &  1.3646e-03  &  0.9846  &  4.8630e-04  &  0.9862  &  1.8793e-03  &  0.9865 
% 3.12500e-04  &  6.7767e-04  &  1.0098  &  2.4133e-04  &  1.0108  &  9.3275e-04  &  1.0107 
% 1.56250e-04  &  3.2923e-04  &  1.0415  &  1.1720e-04  &  1.0421  &  4.5293e-04  &  1.0422  


%%linear BDF2
% ex03_1_Vesicles_data
%    dt 	   & phi_Error_L2 &phi_Order &  u_Error_L2  &  u_Order &  p_Error_L2  &  p_Order 
% 5.00000e-03  &  8.3228e-05  &  0.0000  &  1.3697e-07  &  0.0000  &  5.2095e-06  &  0.0000 
% 2.50000e-03  &  2.1253e-05  &  1.9694  &  3.4589e-08  &  1.9855  &  1.3159e-06  &  1.9852 
% 1.25000e-03  &  5.3732e-06  &  1.9838  &  8.7157e-09  &  1.9886  &  3.3066e-07  &  1.9926 
% 6.25000e-04  &  1.3508e-06  &  1.9919  &  2.1889e-09  &  1.9934  &  8.2867e-08  &  1.9965 
% 3.12500e-04  &  3.3845e-07  &  1.9969  &  5.4821e-10  &  1.9974  &  2.0728e-08  &  1.9992 
% 1.56250e-04  &  8.4476e-08  &  2.0023  &  1.3682e-10  &  2.0025  &  5.1695e-09  &  2.0035 
% 7.81250e-05  &  2.0873e-08  &  2.0169  &  3.3811e-11  &  2.0167  &  1.2768e-09  &  2.0175 
% 3.90625e-05  &  4.9589e-09  &  2.0736  &  8.0375e-12  &  2.0727  &  3.0349e-10  &  2.0728 


% ex03_2_Vesicles_data
% Q=1.0000e+00,eta2=2.5314e+01,epsilon=0.147,t=0.00199/0.0020, dt=3.91e-06, Nx=64, Ny=64, timeElapsed=7.815106
%    dt 	   & phi_Error_L2 &phi_Order &  u_Error_L2  &  u_Order &  p_Error_L2  &  p_Order 
% 5.00000e-04  &  6.3256e-05  &  0.0000  &  6.6450e-04  &  0.0000  &  7.3069e-03  &  0.0000 
% 2.50000e-04  &  1.5545e-05  &  2.0248  &  3.2507e-04  &  1.0315  &  3.5617e-03  &  1.0367 
% 1.25000e-04  &  3.8132e-06  &  2.0274  &  1.6459e-04  &  0.9819  &  4.2453e-03  &  -0.2533 
% 6.25000e-05  &  9.7485e-07  &  1.9677  &  8.2891e-05  &  0.9896  &  3.7791e-03  &  0.1678 
% 3.12500e-05  &  2.6695e-07  &  1.8686  &  4.1540e-05  &  0.9967  &  3.4479e-03  &  0.1323 
% 1.56250e-05  &  8.4688e-08  &  1.6563  &  2.0799e-05  &  0.9980  &  2.3784e-03  &  0.5357 
% 7.81250e-06  &  3.2990e-08  &  1.3601  &  1.0388e-05  &  1.0016  &  2.5502e-03  &  -0.1006 
% 3.90625e-06  &  1.4921e-08  &  1.1446  &  5.1887e-06  &  1.0015  &  1.5827e-03  &  0.6882 

%% nonlinear_1st
% ex03_1_Vesicles_data
% Q=1.0000e+00,eta_1=1.0000e+00,eta_2=-3.3247e+01,epsilon=0.147,t=0.20000/0.2000, dt=3.91e-05, Nx=64, Ny=64, timeElapsed=104.618760
%    dt 	   & phi_Error_L2 &phi_Order &  u_Error_L2  &  u_Order &  p_Error_L2  &  p_Order 
% 5.00000e-03  &  1.8942e-03  &  0.0000  &  5.1315e-06  &  0.0000  &  5.4602e-04  &  0.0000 
% 2.50000e-03  &  1.0109e-03  &  0.9059  &  2.7831e-06  &  0.8827  &  2.9207e-04  &  0.9026 
% 1.25000e-03  &  5.2151e-04  &  0.9549  &  1.4493e-06  &  0.9413  &  1.5052e-04  &  0.9564 
% 6.25000e-04  &  2.6340e-04  &  0.9855  &  7.3570e-07  &  0.9782  &  7.6066e-05  &  0.9846 
% 3.12500e-04  &  1.3075e-04  &  1.0104  &  3.6618e-07  &  1.0066  &  3.7772e-05  &  1.0099 
% 1.56250e-04  &  6.3510e-05  &  1.0418  &  1.7810e-07  &  1.0399  &  1.8350e-05  &  1.0416 
% 7.81250e-05  &  2.9653e-05  &  1.0988  &  8.3210e-08  &  1.0978  &  8.5680e-06  &  1.0987 
% 3.90625e-05  &  1.2665e-05  &  1.2273  &  3.5552e-08  &  1.2269  &  3.6596e-06  &  1.2273 

% ex03_2_Vesicles_data(64)
% Q=1.0000e+00,eta_1=1.0000e+00,eta_2=2.5347e+01,epsilon=0.147,t=0.00199/0.0020, dt=3.91e-06, Nx=64, Ny=64, timeElapsed=11.903718
%    dt 	   & phi_Error_L2 &phi_Order &  u_Error_L2  &  u_Order &  p_Error_L2  &  p_Order 
% 5.00000e-04  &  1.8199e-04  &  0.0000  &  9.7015e-04  &  0.0000  &  1.2942e-02  &  0.0000 
% 2.50000e-04  &  9.7829e-05  &  0.8955  &  5.1107e-04  &  0.9247  &  7.1968e-03  &  0.8466 
% 1.25000e-04  &  5.0775e-05  &  0.9461  &  2.6298e-04  &  0.9586  &  3.8258e-03  &  0.9116 
% 6.25000e-05  &  2.5867e-05  &  0.9730  &  1.3350e-04  &  0.9782  &  1.9774e-03  &  0.9521 
% 3.12500e-05  &  1.3053e-05  &  0.9867  &  6.7263e-05  &  0.9889  &  1.0059e-03  &  0.9752 
% 1.56250e-05  &  6.5553e-06  &  0.9937  &  3.3755e-05  &  0.9947  &  5.0726e-04  &  0.9876 
% 7.81250e-06  &  3.2832e-06  &  0.9975  &  1.6901e-05  &  0.9980  &  2.5461e-04  &  0.9944 
% 3.90625e-06  &  1.6414e-06  &  1.0002  &  8.4481e-06  &  1.0004  &  1.2743e-04  &  0.9986 



%% nonlinear_bdf2
% ex03_1_Vesicles_data
% Q=1.0000e+00,eta_1=1.0000e+00,eta_2=-3.3247e+01,epsilon=0.147,t=0.20000/0.2000, dt=3.91e-05, Nx=64, Ny=64, timeElapsed=127.451341
%    dt 	    & phi_Error_L2  &phi_Order &  u_Error_L2  &  u_Order &  p_Error_L2  &  p_Order 
% 5.00000e-03  &  8.3226e-05  &  0.0000  &  1.3696e-07  &  0.0000  &  5.2074e-06  &  0.0000 
% 2.50000e-03  &  2.1253e-05  &  1.9694  &  3.4585e-08  &  1.9855  &  1.3153e-06  &  1.9852 
% 1.25000e-03  &  5.3729e-06  &  1.9839  &  8.7156e-09  &  1.9885  &  3.3067e-07  &  1.9919 
% 6.25000e-04  &  1.3508e-06  &  1.9919  &  2.1888e-09  &  1.9934  &  8.2870e-08  &  1.9965 
% 3.12500e-04  &  3.3843e-07  &  1.9969  &  5.4820e-10  &  1.9974  &  2.0729e-08  &  1.9992 
% 1.56250e-04  &  8.4472e-08  &  2.0023  &  1.3682e-10  &  2.0025  &  5.1701e-09  &  2.0034 
% 7.81250e-05  &  2.0872e-08  &  2.0169  &  3.3810e-11  &  2.0167  &  1.2774e-09  &  2.0169 
% 3.90625e-05  &  4.9587e-09  &  2.0736  &  8.0374e-12  &  2.0727  &  3.0373e-10  &  2.0724 

% S = 0
% Q=1.0000e+00,eta_1=1.0000e+00,eta_2=-3.3247e+01,epsilon=0.147,t=0.20000/0.2000, dt=3.91e-05, Nx=64, Ny=64, timeElapsed=112.778240
%    dt 	   & phi_Error_L2 &phi_Order &  u_Error_L2  &  u_Order &  p_Error_L2  &  p_Order 
% 5.00000e-03  &  8.0055e-06  &  0.0000  &  1.5801e-08  &  0.0000  &  2.2220e-06  &  0.0000 
% 2.50000e-03  &  2.0088e-06  &  1.9946  &  3.9463e-09  &  2.0015  &  5.5691e-07  &  1.9963 
% 1.25000e-03  &  5.0290e-07  &  1.9980  &  9.8742e-10  &  1.9988  &  1.3952e-07  &  1.9970 
% 6.25000e-04  &  1.2582e-07  &  1.9989  &  2.4680e-10  &  2.0003  &  3.4890e-08  &  1.9996 
% 3.12500e-04  &  3.1445e-08  &  2.0005  &  6.1659e-11  &  2.0010  &  8.7181e-09  &  2.0007 
% 1.56250e-04  &  7.8382e-09  &  2.0042  &  1.5373e-11  &  2.0039  &  2.1732e-09  &  2.0042 
% 7.81250e-05  &  1.9351e-09  &  2.0181  &  3.8012e-12  &  2.0159  &  5.3683e-10  &  2.0173 
% 3.90625e-05  &  4.5917e-10  &  2.0753  &  9.0733e-13  &  2.0668  &  1.2769e-10  &  2.0718 

% ex03_2_Vesicles_data(64)
% Q=1.0000e+00,eta_1=1.0000e+00,eta_2=2.5314e+01,epsilon=0.147,t=0.00199/0.0020, dt=3.91e-06, Nx=64, Ny=64, timeElapsed=11.187836
%    dt 	   & phi_Error_L2 &phi_Order &  u_Error_L2  &  u_Order &  p_Error_L2  &  p_Order 
% 5.00000e-04  &  6.3207e-05  &  0.0000  &  6.6451e-04  &  0.0000  &  7.2833e-03  &  0.0000 
% 2.50000e-04  &  1.5533e-05  &  2.0248  &  3.2507e-04  &  1.0315  &  3.5615e-03  &  1.0321 
% 1.25000e-04  &  3.8110e-06  &  2.0271  &  1.6459e-04  &  0.9819  &  4.2453e-03  &  -0.2534 
% 6.25000e-05  &  9.7474e-07  &  1.9671  &  8.2891e-05  &  0.9896  &  3.7791e-03  &  0.1678 
% 3.12500e-05  &  2.6710e-07  &  1.8676  &  4.1540e-05  &  0.9967  &  3.4479e-03  &  0.1323 
% 1.56250e-05  &  8.4786e-08  &  1.6555  &  2.0799e-05  &  0.9980  &  2.3784e-03  &  0.5357 
% 7.81250e-06  &  3.3028e-08  &  1.3602  &  1.0388e-05  &  1.0016  &  2.5502e-03  &  -0.1006 
% 3.90625e-06  &  1.4933e-08  &  1.1452  &  5.1887e-06  &  1.0015  &  1.5828e-03  &  0.6882 
