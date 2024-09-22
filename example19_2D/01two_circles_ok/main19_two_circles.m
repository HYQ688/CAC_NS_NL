close all;
clear; clc;

% add path
addpath('../','-begin');

% Space: Domain and N
% d = 0.54*pi;
domain.left   =   0 ;
domain.right  =2*pi ;
domain.bottom = -pi ;
domain.top    =  pi ;

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;

N = 64;
Nx = N;
Ny = N;

scheme1 = 'linear';
% scheme1 = 'nonlinear';

% scheme2 = '_1st';   % First-order scheme
% scheme2 = '_2cn';   % Second-order CN scheme
scheme2 = '_bdf2';  % Second-order BDF scheme

scheme = [scheme1 , scheme2];

% Parameters
para.epsilon = 6*pi/128;
para.gamma = 2;
para.Re = 1;
para.lambda = 0.01;
para.C0 = 100;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;

para.name = '_CAC_Vesicles_data_LM_two_circles';
para.name = [scheme,para.name];

% Time: dt T
% T = 0.2;
% t0 = 0;
% tsave = 0.02*T;
% 
% dt_array = 1./2.^(0:4); 
% dt_ref = 1e-6;

T = 8;
t0 = 0;
tsave = 0.02*T;

dt_array = 0.01./2.^(0:4); 
dt_ref = 1e-4;

maxIt = length(dt_array);

%% option
option.plotflag  = 0;
option.printflag = 1;
option.vtkflag  = 0;
option.saveflag  = 1;
option.savefinal  = 0;
option.energyflag = 0;
option.tol = 1e-14;
option.tolit = 1e-11;
option.maxit = 2000;

pde = ex02_Vesicles_data(para);

if 1 == strcmp(scheme,'linear_1st')
    solver_fun = @CAC_Vesicle_2D_linear_newLM_p_SAV_1st;
elseif 1 == strcmp(scheme,'linear_2cn')
    solver_fun = @CAC_Vesicle_2D_linear_newLM_p_SAV_2cn;
elseif 1 == strcmp(scheme,'linear_bdf2')
    solver_fun = @CAC_Vesicle_with_NS_2D_linear_newLM_p_SAV_bdf2;
elseif 1 == strcmp(scheme,'nonlinear_1st')
    solver_fun = @CAC_Vesicle_2D_newLM_p_SAV_1st;
elseif 1 == strcmp(scheme,'nonlinear_2cn')
    solver_fun = @CAC_Vesicle_2D_newLM_p_SAV_2cn;
elseif 1 == strcmp(scheme,'nonlinear_bdf2')
    solver_fun = @CAC_Vesicle_with_NS_2D_newLM_p_SAV_bdf2;
end

%% Run:
if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
    time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
    solver_fun(pde,domain,Nx,Ny,time,option);
end

