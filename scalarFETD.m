% 2D SCALAR FINITE ELEMENT - TIME DOMAIN METHOD WITH FIRST ORDER ABC
% 
clear all; close all; clc; tic;
%% settings
N = 100; % interpolation size
Hmax=0.2; % maximum element size

% physical quantities
c = 3e8;
lambda = 1;
f = c/lambda;
omega = 2*pi*f;
k0 = 2*pi/lambda;
eps0 = 8.854187817e-12;  
mu0 = 4*pi*1e-7;
eta0 = sqrt(mu0/eps0); % wave impedance

% time step setting
dt = Hmax/c/20;
NumTimeSteps = 10000;

a = lambda*2/3; % radius of circle object
d = lambda*2; % distance between obj and truncation box

% the incident angle measured with respect to the positive x-axis
phi_in = 1*pi/4; 

% choose what to plot and how to plot
plot_logistics = true; % if to plot FEM matrices and meshes
plot_inc = false;
use_pdeplot = true;

% (incident field) (not used)
Ez_inc_amp = 1;
Ez_inc_dir = [cos(phi_in),sin(phi_in)];
Jz_inc_amp = 1; tau_p = 3/omega;

%% macros and lambda expressions
Nel2D = @(x,y,ael,bel,cel,Delta_e) 1/(2*Delta_e)*(ael+bel*x+cel*y); 
NablaNel2D = @(bel,cel,Delta_e) 1/(2*Delta_e).*[bel,cel];

GaussPointsLine5 = [-0.9061798459,-0.5384693101,0,0.5384693101,0.9061798459];
GaussWeightsLine5 = [0.2369268851,0.4786286705,0.5688888889,0.4786286705,0.2369268851];
GaussQuadLine5 = @(GaussSamples) GaussSamples*GaussWeightsLine5';

GaussWeightsTri3 = [1/3,1/3,1/3];
GaussQuadTri3 = @(GaussSamples) GaussSamples*GaussWeightsTri3';

Rotate90Mat = [cosd(90) -sind(90); sind(90) cosd(90)];
NextLocalNodeof = [2,3,1]; % map to the next
PreviousLocalNodeof = [3,1,2]; % map to the previous
Ezinc = @(x,y,t) Ez_inc_amp*exp(-1j*k0*dot([x,y],Ez_inc_dir))*exp(1j*omega*t);
dEzinc_dt = @(x,y,t) 1j*omega*Ez_inc_amp*exp(-1j*k0*dot([x,y],Ez_inc_dir))*exp(1j*omega*t);
Jzinc = @(t) Jz_inc_amp*(1-exp(-(t)/tau_p)).*sin(omega*t);
dJzinc_dt = @(t) omega*Jz_inc_amp*(1-exp(-(t)/tau_p)).*cos(omega*t)+Jz_inc_amp*(exp(-(t)/tau_p)/tau_p).*sin(omega*t);

%% Visualization of the incident field (plane wave)
X=linspace(-a-d,a+d,N);
Y=linspace(-a-d,a+d,N);
EzIncGrid=zeros(N,N); 
if plot_inc
    for n=1:NumTimeSteps
        t = n*dt;
        for i=1:N
            for j=1:N
                x = X(i);
                y = Y(j);
                [inPEC,onPEC]=incircle([x,y],[0,0],a);

                if ~inPEC && ~onPEC
                    EzIncGrid(i,j) = Ezinc(x,y,t);
                end
            end
        end
        figure(1)
        imagesc(X,Y,abs(real(EzIncGrid)));
        colormap summer
        axis image
        xlabel('Y (m)')
        ylabel('X (m)')
        title(['$E_z^{inc}$', ' at ', 'time = ',num2str(round(n*dt,9)),' s',],'interpreter','latex');
        set(gca,'fontsize',18);
        drawnow;
        pause(0.1);	
    end
end
%% Geometries and materials
model = createpde;
Circle = [1 0 0 a]';
Rect = [3 4 -a-d -a-d a+d a+d -a-d a+d a+d -a-d]';
Circle = [Circle;zeros(length(Rect) - length(Circle),1)];
gd = [Circle, Rect];
ns = char('Circle','Rect');
ns = ns';
sf = 'Rect-Circle';
g= decsg(gd,sf,ns);
geometryFromEdges(model,g);
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',Hmax);

if plot_logistics
    figure; hold on;
    pdegplot(model,'FaceLabels','on','EdgeLabels','on')
    xlabel('x');
    xlabel('y');
    hold off;

    figure; hold on;
    pdemesh(mesh,'ElementLabels','off')
    xlabel('x');
    xlabel('y');
    set(gca,'fontsize',24);
    hold off;
end

Num_Nodes = size(mesh.Nodes,2);
Num_Elements = size(mesh.Elements,2);

eps = ones(1,Num_Elements); % relative permittivity of the elements
mu = ones(1,Num_Elements); % relative permeability of the elements
sigma = zeros(1,Num_Elements); % conductivity

%% Identify the source and bounaries
NodesOnPEC = findNodes(mesh,'region','Edge',[5, 6, 7, 8]);
NodesOnBoundary = findNodes(mesh,'region','Edge',[1, 2, 3, 4]);
NodeNearSource = findNodes(mesh,'nearest',[-1;-1]);
    
%% Assembly of the spatial matrix system.
T = zeros(Num_Nodes,Num_Nodes);
R = zeros(Num_Nodes,Num_Nodes);
S = zeros(Num_Nodes,Num_Nodes);

xe = zeros(1,3);
ye = zeros(1,3);
ae = zeros(1,3);
be = zeros(1,3);
ce = zeros(1,3);
    
for e=1:Num_Elements
    % initialization of the local element matrices
    Te = zeros(3,3);
    Re = zeros(3,3);
    Se = zeros(3,3);
    
    ne(1:3) = mesh.Elements(1:3,e); % the global node idx of node 1 in element e
    for i=1:3 % the global coordinate of node 1 in element e
        where_am_I = mesh.Nodes(:,ne(i)); 
        xe(i) = where_am_I(1);
        ye(i) = where_am_I(2);
    end
    
    eps_e = eps0*eps(e);
    mu_e = mu0*mu(e);
    sigma_e = 1*sigma(e);
    eta_e = eta0*sqrt(mu(e)/eps(e));
    
    ae(1) = xe(2)*ye(3)-xe(3)*ye(2);
    ae(2) = xe(3)*ye(1)-xe(1)*ye(3);
    ae(3) = xe(1)*ye(2)-xe(2)*ye(1);
    be(1) = ye(2) - ye(3);% y 2-3
    be(2) = ye(3) - ye(1);% y 3-1
    be(3) = ye(1) - ye(2);% y 1-2
    ce(1) = xe(3) - xe(2);% x 3-2
    ce(2) = xe(1) - xe(3);% x 1-3
    ce(3) = xe(2) - xe(1);% x 2-1
    Delta_e = (be(1)*ce(2)-be(2)*ce(1))/2;
    
    for i=1:3
        for j=1:3
            NablaNei=NablaNel2D(be(i),ce(i),Delta_e);
            NablaNej=NablaNel2D(be(j),ce(j),Delta_e);
            Se(i,j)=Se(i,j)+1/mu_e*dot(NablaNei,NablaNej)*Delta_e;
            Te(i,j)=Te(i,j)+eps_e*Delta_e/12*(1+(i==j));
            Re(i,j)=Re(i,j)+sigma_e*Delta_e/12*(1+(i==j));
        end
    end
    
    for lnode = 1:3 % local node
        nlnode = NextLocalNodeof(lnode); % next local node
        plnode = PreviousLocalNodeof(lnode); % previous local node

        % if either edge 1-2, 2-3, 3-1 is on the outer boundary so that the edge is not shared by any neighboring two elements (not a common edge so that their contribution will not be canceled)
        if ismember(ne(lnode),NodesOnBoundary) && ismember(ne(nlnode),NodesOnBoundary)
            vecEdge = [xe(nlnode)-xe(lnode),ye(nlnode)-ye(lnode)];
            edgeLength = norm(vecEdge);
            for i=[lnode, nlnode]
                for j=[lnode, nlnode]
                    Re(i,j)=Re(i,j)+1/eta_e*edgeLength/6*(1+(i==j));
                end
            end
        end

        if ismember(ne(lnode),NodesOnBoundary) && ismember(ne(plnode),NodesOnBoundary)
            vecEdge = [xe(plnode)-xe(lnode),ye(plnode)-ye(lnode)];
            edgeLength = norm(vecEdge);
            for i=[plnode, lnode]
                for j=[plnode, lnode]
                    Re(i,j)=Re(i,j)+1/eta_e*edgeLength/6*(1+(i==j));
                end
            end
        end
    end
    
    % Add local to global
    for i=1:3
        for j=1:3
            T(ne(i),ne(j)) = T(ne(i),ne(j))+Te(i,j); 
            R(ne(i),ne(j)) = R(ne(i),ne(j))+Re(i,j); 
            S(ne(i),ne(j)) = S(ne(i),ne(j))+Se(i,j); 
        end
    end
end

%% Time-stepping and Solving Linear System
K = 1/dt^2*T+1/2/dt*R;
EzTot_coefficients1 = zeros(Num_Nodes,1);
EzTot_coefficients2 = zeros(Num_Nodes,1);
EzTot_coefficients3 = zeros(Num_Nodes,1);

for n=1:NumTimeSteps
    t = n*dt;
    l = zeros(Num_Nodes,1);
    for e=1:Num_Elements
        le = zeros(1,3);
        
        ne(1:3) = mesh.Elements(1:3,e); % the global node idx of node 1 in element e
        for i=1:3 % the global coordinate of node 1 in element e
            where_am_I = mesh.Nodes(:,ne(i)); 
            xe(i) = where_am_I(1);
            ye(i) = where_am_I(2);
        end

        eps_e = eps0*eps(e);
        mu_e = mu0*mu(e);
        sigma_e = 1*sigma(e);
        eta_e = eta0*sqrt(mu(e)/eps(e));

        ae(1) = xe(2)*ye(3)-xe(3)*ye(2);
        ae(2) = xe(3)*ye(1)-xe(1)*ye(3);
        ae(3) = xe(1)*ye(2)-xe(2)*ye(1);
        be(1) = ye(2) - ye(3);% y 2-3
        be(2) = ye(3) - ye(1);% y 3-1
        be(3) = ye(1) - ye(2);% y 1-2
        ce(1) = xe(3) - xe(2);% x 3-2
        ce(2) = xe(1) - xe(3);% x 1-3
        ce(3) = xe(2) - xe(1);% x 2-1
        Delta_e = (be(1)*ce(2)-be(2)*ce(1))/2;
        
        % evaluate the RHS of the 2nd-order updating equations
        for lnode = 1:3 % local node
            nlnode = NextLocalNodeof(lnode); % next local node
            plnode = PreviousLocalNodeof(lnode); % previous local node
            
            if ne(lnode)==NodeNearSource
                GaussSamples = zeros(1,3);
                for ngp=1:3
                    % 1. tranform the quadrature points from natural coordinates to
                    % cartesian coordinates (3 points Gauss quadrature, 2nd order
                    % precision) 
                    x = (xe(nlnode)+xe(plnode))/2;
                    y = (ye(nlnode)+ye(plnode))/2;
                    % 2. evaluate functions and their dot product at the points
                    Ne_i = Nel2D(x,y,ae(lnode),be(lnode),ce(lnode),Delta_e);
                    % 3. multiply with Gauss weights and scale
                    GaussSamples(ngp) = Ne_i*dJzinc_dt(t);
                end
                le(lnode)=le(lnode)-GaussQuadTri3(GaussSamples)*2*Delta_e;
            end
        end
        for i=1:3
            l(ne(i))=l(ne(i))+le(i);
        end
    end
    b = (2/dt^2*T-S)*EzTot_coefficients2-(1/dt^2*T-1/2/dt*R)*EzTot_coefficients1+l;
    for i=1:length(NodesOnPEC) % Dirichlet on PEC
        b(NodesOnPEC(i))=0;
        K(NodesOnPEC(i),NodesOnPEC(i))=1;
        for j=1:Num_Nodes
            if j==NodesOnPEC(i), continue; end
            K(NodesOnPEC(i),j)=0;
            K(j,NodesOnPEC(i))=0;
        end
    end
    
    
    EzTot_coefficients3 = K\b;
    EzTot_coefficients1 = EzTot_coefficients2;
    EzTot_coefficients2 = EzTot_coefficients3;
    
    if use_pdeplot
        if mod(n,10)==0
            figure(5);
            pdeplot(model, 'XYData', abs(real(EzTot_coefficients3)));
            colormap summer;colorbar;axis image;xlabel('Y (m)');ylabel('X (m)');set(gca,'fontsize',18);
            title(['$E_z^{tot}$', ' at ', 'time = ',num2str(round(n*dt,9)),' s',],'interpreter','latex');
            drawnow; pause(0.2);
        end
    else
        if mod(n,50)==0
            EzGrid = zeros(N,N);
            GridInterpolated = zeros(N,N);
            for e=1:Num_Elements

                ne(1:3) = mesh.Elements(1:3,e); % the global node idx of node 1 in element e
                for i=1:3 % the global coordinate of node 1 in element e
                    where_am_I = mesh.Nodes(:,ne(i)); 
                    xe(i) = where_am_I(1);
                    ye(i) = where_am_I(2);
                end

                ae(1) = xe(2)*ye(3)-xe(3)*ye(2);
                ae(2) = xe(3)*ye(1)-xe(1)*ye(3);
                ae(3) = xe(1)*ye(2)-xe(2)*ye(1);
                be(1) = ye(2) - ye(3);% y 2-3
                be(2) = ye(3) - ye(1);% y 3-1
                be(3) = ye(1) - ye(2);% y 1-2
                ce(1) = xe(3) - xe(2);% x 3-2
                ce(2) = xe(1) - xe(3);% x 1-3
                ce(3) = xe(2) - xe(1);% x 2-1
                Delta_e = (be(1)*ce(2)-be(2)*ce(1))/2;
                for q=1:N
                    for p=1:N
                        if ~GridInterpolated(q,p)
                            [in,on]=inpolygon(X(p),Y(q),xe,ye);
                            if (in || on) %indicating if inside the element or on the edge of the element
                                for i=1:3
                                    EzGrid(q,p)=EzGrid(q,p)+Nel2D(X(p),Y(q),ae(i),be(i),ce(i),Delta_e)*EzTot_coefficients3(ne(i));
                                end
                                GridInterpolated(q,p) = 1;
                            end
                        end
                    end
                end
            end
            figure(5);
            imagesc(X,Y,abs(real(EzGrid)));
            colormap summer;colorbar;axis image;xlabel('Y (m)');ylabel('X (m)');set(gca,'fontsize',18);
            title(['$E_z^{tot}$', ' at ', 'time = ',num2str(round(n*dt,9)),' s',],'interpreter','latex');
            drawnow; 
            set(gcf,'Position',[0,0,1024,1024]); 
            saveas(gca, ['step',num2str(n),'Ez','.png']);
            pause(0.2);
        end
    end
end
toc;
%% functions
function [in,on] = incircle(coord,center,radius)
    in = norm(coord-center)<radius;
    on = norm(coord-center)==radius;
end