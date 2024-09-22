function pde = ex03_1_Vesicles_data(para)
if nargin == 0
    epsilon = 1;
    gamma = 1;
    lambda  = 1;
    C0 = 0;
    Re = 1;
    S1 = 0;
    S2 = 0;
    S3 = 0;
    name = '';
else
    if ~isstruct(para)
        exit('we need a struct data');
    end
    if ~isfield(para,'epsilon') || isempty(para.epsilon)
        epsilon = 1;
    else
        epsilon = para.epsilon;
    end
    if ~isfield(para,'gamma') || isempty(para.gamma)
        gamma = 1;
    else
        gamma = para.gamma;
    end
    if ~isfield(para,'Re') || isempty(para.Re)
        Re = 1;
    else
        Re = para.Re;
    end 
    if ~isfield(para,'lambda') || isempty(para.lambda)
        lambda = 1;
    else
        lambda = para.lambda;
    end
    if ~isfield(para,'C0') || isempty(para.C0)
        C0 = 0;
    else
        C0 = para.C0;
    end
    if ~isfield(para,'S1') || isempty(para.S1)
        S1 = 0;
    else
        S1 = para.S1;
    end
    if ~isfield(para,'S2') || isempty(para.S2)
        S2 = 0;
    else
        S2 = para.S2;
    end
    if ~isfield(para,'S3') || isempty(para.S3)
        S3 = 0;
    else
        S3 = para.S3;
    end
    if ~isfield(para,'name') || isempty(para.name)
        name = 'ex03_1_Vesicles_data';
    else
        name = para.name;
    end
end

pde = struct('epsilon',epsilon, ...
    'gamma',gamma, ...
    'lambda',lambda, ...
    'Re',Re, ...
    'C0',C0, ...
    'S1',S1, ...
    'S2',S2, ...
    'S3',S3, ...
    'initphi',@initphi, ...
    'initu',@initu,...
    'initv',@initv,...
    'initp',@initp,...
    'exactphi',@exactphi, ...
    'exactu',@exactu,...
    'exactv',@exactv,...
    'exactp',@exactp,...
    'rhs_u',@rhs_u,...
    'rhs_v',@rhs_v,...
    'name',name);

    function z = initphi(x,y)
        z = exactphi(x,y,0);
    end

    function z = exactphi(x,y,t)
%         z = (sin(2*x).*sin(2*y)/4+0.48).*(1-sin(t).^2/2);
        z = exp(-t).*sin(2.*x).*cos(2.*y)./4;
    end

    function z = initu(x,y)
        z = exactu(x,y,0);
    end

    function z = initv(x,y)
        z = exactv(x,y,0);
    end

    function z = initp(x,y)
        z = exactp(x,y,0);
    end

    function z = exactu(x,y,t)
        z = sin(2*y).*sin(x).^2.*sin(t);
    end

    function z = exactv(x,y,t)
        z = -sin(2*x).*sin(y).^2.*sin(t); 
    end

    function z = exactp(x,y,t)
        z = cos(x).*sin(t).*sin(y);
    end

%     function rhs = rhs_u(x,y,t)
%         rhs = sin(2.*y).*cos(t).*sin(x).^2 - sin(t).*sin(x).*sin(y) + 2.*sin(2.*y).^2.*cos(x).*sin(t).^2.*sin(x).^3 - (2.*sin(2.*y).*cos(x).^2.*sin(t))./Re + (6.*sin(2.*y).*sin(t).*sin(x).^2)./Re + 2.*lambda.*cos(2.*x).*cos(2.*y).^2.*sin(2.*x).*exp(-2.*t) - 4.*cos(2.*y).*cos(x).*sin(t).^2.*sin(x).^3.*sin(y).^2 - (3.*lambda.*cos(2.*x).*cos(2.*y).^4.*sin(2.*x).^3.*exp(-4.*t))./(8.*epsilon) + (3.*lambda.*cos(2.*x).^3.*cos(2.*y).^4.*sin(2.*x).*exp(-4.*t))./(16.*epsilon) + (lambda.*cos(2.*x).*cos(2.*y).^4.*sin(2.*x).^3.*exp(-4.*t))./(32.*epsilon.^3) - (3.*lambda.*cos(2.*x).*cos(2.*y).^6.*sin(2.*x).^5.*exp(-6.*t))./(2048.*epsilon.^3) - 8.*epsilon.*lambda.*cos(2.*x).*cos(2.*y).^2.*sin(2.*x).*exp(-2.*t) - (lambda.*cos(2.*x).*cos(2.*y).^2.*sin(2.*x).*exp(-2.*t))./(8.*epsilon.^3) + (3.*lambda.*cos(2.*x).*cos(2.*y).^2.*sin(2.*x).^3.*sin(2.*y).^2.*exp(-4.*t))./(16.*epsilon);
%     end
% 
%     function rhs = rhs_v(x,y,t)
%         rhs = cos(x).*cos(y).*sin(t) - sin(2.*x).*cos(t).*sin(y).^2 + 2.*sin(2.*x).^2.*cos(y).*sin(t).^2.*sin(y).^3 + (2.*sin(2.*x).*cos(y).^2.*sin(t))./Re - (6.*sin(2.*x).*sin(t).*sin(y).^2)./Re - 2.*lambda.*cos(2.*y).*sin(2.*x).^2.*sin(2.*y).*exp(-2.*t) - 4.*cos(2.*x).*cos(y).*sin(t).^2.*sin(x).^2.*sin(y).^3 - (3.*lambda.*cos(2.*y).*sin(2.*x).^4.*sin(2.*y).^3.*exp(-4.*t))./(16.*epsilon) + (3.*lambda.*cos(2.*y).^3.*sin(2.*x).^4.*sin(2.*y).*exp(-4.*t))./(8.*epsilon) - (lambda.*cos(2.*y).^3.*sin(2.*x).^4.*sin(2.*y).*exp(-4.*t))./(32.*epsilon.^3) + (3.*lambda.*cos(2.*y).^5.*sin(2.*x).^6.*sin(2.*y).*exp(-6.*t))./(2048.*epsilon.^3) + 8.*epsilon.*lambda.*cos(2.*y).*sin(2.*x).^2.*sin(2.*y).*exp(-2.*t) + (lambda.*cos(2.*y).*sin(2.*x).^2.*sin(2.*y).*exp(-2.*t))./(8.*epsilon.^3) - (3.*lambda.*cos(2.*x).^2.*cos(2.*y).^3.*sin(2.*x).^2.*sin(2.*y).*exp(-4.*t))./(16.*epsilon);
%     end


end