% 2d scalar (nodal) FEM formulation for scattering problem using
% PML. The finite elements are triangles with linear basis
% fucntions. Objects are PECs. Ez^sca is the primary unknown
close all; clear all; clc;
%% parameters and settings
R = 5;
c = 3e8; % speed of light
lambda=c/3e8;
% angular frequency or wave number, we can use either one of them to
% formulate the wave propagation problem:
% 1. if we use angular frequency, the material properties should be
% expressed in absolute form
% 2. if we use wave number, the material properties can be relative ones
k0 = 2*pi/lambda;
f=c/lambda; % lambda = c/f
omega =2*pi*f;
E0 = 1;
phi_in = 1*pi/4;
k_dir = [cos(phi_in),sin(phi_in)]; % wave vector
plot_logistics = true;
a = lambda; % radius of circle object
d = lambda; % distance between obj and pml

eps0 = 8.854187817e-12;  
mu0 = 4*pi*1e-7;
eta0 = sqrt(mu0/eps0); % wave impedance

%% Macros
Ez_inc = @(x,y) E0*exp(-1j*k0*([x,y]*k_dir'));
Nel2D = @(x,y,ael,bel,cel,Delta_e) ael+bel*x+cel*y; 
NablaNel2D = @(bel,cel) [bel,cel];
%% PML 
m = 0; % mth-order polynomial
L = lambda/2; % PML thickness
pml_attenuation = 1e-6; % specified reflection coefficient
sigma_max = -(m+1)/(2*L*2*pi)*log(pml_attenuation); % maximum conductivity inside the PML
a1=complex(1,-sigma_max);

%% Geometries
C1 = [1 0 0 a]';
rect2 = [3 4 -a-d -a-d a+d a+d -a-d a+d a+d -a-d]';
rect3 = [3 4 -a-d-L -a-d-L -a-d -a-d -a-d a+d a+d -a-d]';
rect4 = [3 4 -a-d-L -a-d-L -a-d -a-d a+d a+d+L a+d+L a+d]';
rect5 = [3 4 -a-d-L -a-d-L -a-d -a-d -a-d-L -a-d -a-d -a-d-L]';
rect6 = [3 4 a+d a+d a+d+L a+d+L -a-d a+d a+d -a-d]';
rect7 = [3 4 a+d a+d a+d+L a+d+L a+d a+d+L a+d+L a+d]';
rect8 = [3 4 a+d a+d a+d+L a+d+L -a-d-L -a-d -a-d -a-d-L]';
rect9 = [3 4 -a-d -a-d a+d a+d a+d a+d+L a+d+L a+d]';
rect10 = [3 4 -a-d -a-d a+d a+d -a-L-d -a-d -a-d -a-d-L]';
C1 = [C1;zeros(length(rect2) - length(C1),1)];
gd = [C1, rect2, rect3, rect4, rect5, rect6, rect7, rect8, rect9, rect10];
ns = char('C1','rect2','rect3','rect4','rect5','rect6','rect7','rect8', 'rect9', 'rect10');
ns = ns';
sf = '(rect2+rect3+rect4+rect5+rect6+rect7+rect8+rect9+rect10)-C1';
dl= decsg(gd,sf,ns);
[p,e,t] = initmesh(dl); % [node, edge, element]
for i=1:1 % number of mesh refinements
    [p,e,t] = refinemesh(dl,p,e,t); 
end

if plot_logistics == 1
    figure; hold on;
    pdegplot(dl,'FaceLabels','on','EdgeLabels','on')
    xlabel('x');
    xlabel('y');
    hold off;

    figure; hold on;
    pdemesh(p,e,t);
    xlabel('x');
    xlabel('y');
    set(gca,'fontsize',24);
    hold off;
end

Num_Nodes = size(p,2); % number of nodes
Num_Edges = size(e,2);  % number of edges
Num_Elements = size(t,2); % number of elements�

eps = ones(1,Num_Elements); % relative permittivity of the elements
mu = ones(1,Num_Elements); % relative permeability of the elements
sigma = zeros(1,Num_Elements); % conductivity

Ez_tot=zeros(Num_Nodes,1); % solution of Kx=b, solved coefficient vector
NodesOnPEC = zeros(Num_Nodes,1); % flag for nodes with unknown coefficients

%         %Find the nodes associated with edges 25, 26, 27 and 28 (the boundary of circle object).
%         Nodes_onObj_Circle = findNodes(mesh,'region','Edge',[25, 26, 27, 28]);
%         Nodes_PEC = Nodes_onObj_Circle;

% excluding all those on the perfectly conducting surface (on PEC boundary of the circle)
for nEdge = 1:Num_Edges
    % e(6,k) is the subdomain number on the left side of the edge.
    % The direction along the edge is given by increasing parameter
    % values. The subdomain 0 is the exterior of the geometry.  
    % e(7,k) is the subdomain number on the right side of the edge. 
    if (e(6,nEdge)==0 || e(7,nEdge)==0) && (e(6,nEdge)==5 || e(7,nEdge)==5)
        NodesOnPEC(e(1,nEdge))=1;
        NodesOnPEC(e(2,nEdge))=1;  
    end
end

K = zeros(Num_Nodes,Num_Nodes); % global LHS matrix
b = zeros(Num_Nodes,1);   % global RHS vector

for ie = 1:Num_Elements          % loop through all elements
    
    Te = zeros(3,3);
    Se = zeros(3,3);
    Ke = zeros(3,3);
    ae=zeros(1,3);
    be=zeros(1,3);
    ce=zeros(1,3);
    
    ne = t(1:3,ie);         % global node number
    sdm = t(4,ie);               % which subdomain(face) the element resides
    xe = p(1,ne(1:3));
    ye = p(2,ne(1:3));
    De = det([1 xe(1) ye(1);1 xe(2) ye(2);1 xe(3) ye(3)]); 
    Ae = abs(De/2);          % Element Area
    % different definition of ae, be and ce, thus different definitions for
    % the basis function
    ae(1) = (xe(2)*ye(3)-xe(3)*ye(2))/De;
    ae(2) = (xe(3)*ye(1)-xe(1)*ye(3))/De;
    ae(3) = (xe(1)*ye(2)-xe(2)*ye(1))/De;
    be(1) = (ye(2)-ye(3))/De;
    ce(1) = (xe(3)-xe(2))/De;
    be(2) = (ye(3)-ye(1))/De;
    ce(2) = (xe(1)-xe(3))/De; 
    be(3) = (ye(1)-ye(2))/De;     
    ce(3) = (xe(2)-xe(1))/De;   
    
    % Assemble local matrices
    for i=1:3 
        for j=1:3 
            NablaNei=NablaNel2D(be(i),ce(i));
            NablaNej=NablaNel2D(be(j),ce(j));
            if sdm==5  %
                ezz=eps0;
                mxx=mu0;
                myy=mu0;
            elseif sdm==3 || sdm==8 || sdm==1 || sdm==9  % for the four corner regions, we choose
                ezz= a1^2*eps0;
                mxx=mu0;
                myy=mu0;
            elseif sdm==4 ||  sdm==6 % for the PML perpendicular to the x-axis, we choose
                ezz= a1*eps0;
                mxx=mu0/a1;
                myy=mu0*a1;
            elseif sdm==2 || sdm==7 % for the PML perpendicular to the y-axis, we choose
                ezz= a1*eps0;
                mxx=a1*mu0;
                myy=mu0/a1;
            end
            
            % PML tensor
            % Lambda = diag([a b c]);
            mu_pml = [myy, mxx];
            eps_pml = ezz;
            Se(i,j) = sum(1./mu_pml.*NablaNei.*NablaNej)*Ae; 
            Te(i,j) = eps_pml*Ae/12*(1+(i==j));
            Ke(i,j) = Se(i,j)-omega^2*Te(i,j);
        end
    end
            
    % Add local matrices to global matrix and vector
    for i=1:3 
        for j=1:3 
            % Quote: At the boundary of the scatter,
            % non-homogeneous Dirichlet conditions are applied. 
            % Therefore at each node of the scatter's boundary, 
            % the scattering field will be known and equal to the incident.  
            if (NodesOnPEC(ne(i))==false) % current node is not on PEC boundary
                if (NodesOnPEC(ne(j))==false) % the edge have no points on PEC boundary
                    K(ne(i),ne(j)) = K(ne(i),ne(j)) + Ke(i,j);
                    b(ne(i)) = b(ne(i));% this case has no contribution to the rhs vector
                else % the other point on this edge is on PEC boundary
                    K(ne(i),ne(j)) = K(ne(i),ne(j)); % this case has no contribution to the lhs matrix
                    b(ne(i)) = b(ne(i))-Ke(i,j)*(-Ez_inc(p(1,ne(j)),p(2,ne(j)))); 
                end
            else % current node is right on PEC boundary
                % the solution to the expansion coefficient of the
                % scattered field Ez^sca on PEC surface is known by applying
                % the physical knowledge of the reflection on the PEC
                % boundary 
                K(ne(i),ne(i)) = 1;
                b(ne(i)) = -Ez_inc(p(1,ne(i)),p(2,ne(i)));
            end
        end
    end
end

% solve
Ez_sca=K\b;

for i=1:Num_Nodes % loop through all the nodes (both known and unknown)
    Ez_tot(i)= Ez_inc(p(1,i),p(2,i)) + Ez_sca(i); % Ez^tot = Ez^inc + Ez^sca
end

% plot
figure; hold on;
pdeplot(p, e, t, 'XYData', abs(Ez_tot));
str = sprintf('E_z^{tot} PML Thickness = %1.2f ,obj_radius = %1.2f', L, a);
title(str)
set(gca,'fontsize',24);
colormap summer
axis image
colorbar;
hold off;

if plot_logistics == 1
    figure; hold on;
    pdeplot(p, e, t, 'XYData', abs(Ez_sca));
    str = sprintf('E_z^{sca} PML Thickness = %1.2f ,obj_radius = %1.2f', L, a);
    title(str)
    set(gca,'fontsize',24);
    colormap summer
    axis image
    colorbar;
    hold off;
end