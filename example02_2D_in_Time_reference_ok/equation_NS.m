%% Epitaxy
clear 
clc
syms x y epsilon t Lx Ly beta gamma ...
    Re lambda M

phi = exp(-t).*sin(2.*x).*cos(2.*y)./4;
u =  sin(2*y).*sin(x).^2.*sin(t);
v = -sin(2*x).*sin(y).^2.*sin(t); 
p  = cos(x).*sin(y).*sin(t); 

div_u = simplify(diff(u,x,1)+ diff(v,y,1));

lap_phi = diff(phi,x,2) + diff(phi,y,2);
lap_u = diff(u,x,2) + diff(u,y,2);
lap_v = diff(v,x,2) + diff(v,y,2);

grad_phi_square = diff(phi,x,1).^2 + diff(phi,y,1).^2;

u1 = diff(u,x,1); u2 = diff(u,y,1);
v1 = diff(v,x,1); v2 = diff(v,y,1);
p1 = diff(p,x,1); p2 = diff(p,y,1);
phi1 = diff(phi,x,1); phi2 = diff(phi,y,1);

div_term = diff(phi.^2.*phi1,x,1) + diff(phi.^2.*phi2,y,1);

convective1 = u.*u1 + v.*u2;
convective2 = u.*v1 + v.*v2;

ut = simplify( diff(u,t));
vt = simplify( diff(v,t));

f     = phi.^3 - phi;
f_der = 3*phi.^2 - 1;

phi_t = simplify(diff(phi,t));

mu = lambda.*epsilon.*( diff(lap_phi,x,2) + diff(lap_phi,y,2)  + 6./epsilon.^2.*(grad_phi_square.*phi - div_term)  + 1./epsilon.^4.*f_der.*f + 2./epsilon.*lap_phi);

f1 = simplify( ut  + convective1 + p1 - 1/Re.*lap_u - (mu.*diff(phi,x,1)));
f2 = simplify( vt  + convective2 + p2 - 1/Re.*lap_v - (mu.*diff(phi,y,1)));


f1 = char(f1); f1 = strrep(f1,'*','.*'); f1 = strrep(f1,'/','./');
f1 = strrep(f1,'^','.^')

f2 = char(f2); f2 = strrep(f2,'*','.*'); f2 = strrep(f2,'/','./');
f2 = strrep(f2,'^','.^')

phi 

phi_t = char(phi_t);
phi_t = strrep(phi_t,'*','.*');
phi_t = strrep(phi_t,'/','./');
phi_t = strrep(phi_t,'^','.^')

u
ut = char(ut);
ut = strrep(ut,'*','.*');
ut = strrep(ut,'/','./');
ut = strrep(ut,'^','.^')

v
vt = char(vt);
vt = strrep(vt,'*','.*');
vt = strrep(vt,'/','./');
vt = strrep(vt,'^','.^')

p
p1
p2
