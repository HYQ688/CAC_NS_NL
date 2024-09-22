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

% Parameters
para.C0 = 10000; % SAV
para.epsilon = 6*pi/128;
% para.epsilon = 0.07;
para.M = 1;
% para.beta_m = 1;
% para.M2 = 1e3;
para.name = 'ex19_BDF2_Vesicles_data_LM_six_circles';

% Time: dt T
% T = 0.2;
% t0 = 0;
% tsave = 0.02*T;
% 
% dt_array = 1./2.^(0:4); 
% dt_ref = 1e-6;

T = 2;
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

%% Run:
if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
    time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
    CAC_Vesicle_2D_BDF2_newLM_p_SAV(pde,domain,Nx,Ny,time,option);
end

