% 2D SCALAR DISCONTINUOUS GALERKIN - TIME DOMAIN METHOD (CENTRAL FLUX) FOR
% TM WAVE WITH FIRST ORDER SILVER MULLER ABC AND LINE CURRENT 
clear all; close all; clc; tic; toc;
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
dt = Hmax/c/10;
NumTimeSteps = 10000;

a = lambda*2/3; % radius of circle object
d = lambda*2; % distance between obj and truncation box
X=linspace(-a-d,a+d,N);
Y=linspace(-a-d,a+d,N);

% choose what to plot and how to plot
plot_logistics = true; % if to plot FEM matrices and meshes
use_pdeplot = false; % if to use the function provided by PDE toolbox for visualization 

% single freq electric current Jz
Jz_inc_amp = 1; tau_p = 3/omega;

%% definitions, macros and lambda expressions
% assumption: when looping through the DOFs in a local element, the loop
% index, in it self without special announcement, is the local node id 

i_map = [3,1,2]; % from loop index to local edge id
l_map = [1,2,3]; % from loop index to local start node id
k_map = [2,3,1]; % from loop index to local end node id
NextLocalNodeof = [2,3,1]; % map to the next
PreviousLocalNodeof = [3,1,2]; % map to the previous

GaussWeightsTri3 = [1/3,1/3,1/3];
GaussQuadTri3 = @(GaussSamples) GaussSamples*GaussWeightsTri3';

Nel2D = @(x,y,ael,bel,cel,Delta_e) 1/(2*Delta_e)*(ael+bel*x+cel*y); 
NablaNel2D = @(bel,cel,Delta_e) 1/(2*Delta_e).*[bel,cel];
Rotate90Mat = [cosd(90) -sind(90); sind(90) cosd(90)];

Jzimp = @(t) Jz_inc_amp*(1-exp(-(t)/tau_p)).*sin(omega*t);
dJzimp_dt = @(t) omega*Jz_inc_amp*(1-exp(-(t)/tau_p)).*cos(omega*t)+Jz_inc_amp*(exp(-(t)/tau_p)/tau_p).*sin(omega*t);

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

% time-dependent solution
ez_coeff_global = zeros(Num_Nodes,1);
ez_coeff1 = zeros(Num_Elements,3);
hx_coeff1 = zeros(Num_Elements,3);
hy_coeff1 = zeros(Num_Elements,3);
ez_coeff2 = zeros(Num_Elements,3);
hx_coeff2 = zeros(Num_Elements,3);
hy_coeff2 = zeros(Num_Elements,3);

eps = ones(1,Num_Elements); % relative permittivity of the elements
mu = ones(1,Num_Elements); % relative permeability of the elements
sigma = zeros(1,Num_Elements); % conductivity

%% Edge numbering and its connection to nodes and elements
% the edges in nodal dgtd are the interfaces between elements where fluxes
% are defined
% edge and its direction (start and end) are uniquely defined globally
Num_Edges = 0;
Element2Edge = zeros(Num_Elements,3);
Edge2Node = zeros(3*Num_Elements,2); % [global_edge_id,global_node_id], do not have a specific direction
Element2SignOfEdge = zeros(Num_Elements,3); % direction

for e=1:Num_Elements
    ne(1:3) = mesh.Elements(1:3,e);
    for i=1:3
        start_node = ne(l_map(i));
        end_node = ne(k_map(i));
        edge_id = find( ( Edge2Node(:,1)==start_node & Edge2Node(:,2)==end_node ) | ( Edge2Node(:,2)==start_node & Edge2Node(:,1)==end_node ), 1);
        
        if isempty(edge_id) % new edge
            % first add it to edge list
            Num_Edges = Num_Edges + 1;
            Edge2Node(Num_Edges,1) = start_node;
            Edge2Node(Num_Edges,2) = end_node;
            % second add the edge to the element-edge list
            Element2Edge(e,i_map(i)) = Num_Edges;
            Element2SignOfEdge(e,i_map(i)) = 1;
        elseif Edge2Node(edge_id,2)==start_node && Edge2Node(edge_id,1)==end_node % shared edge with opposite direction
            Element2Edge(e,i_map(i)) = edge_id;
            Element2SignOfEdge(e,i_map(i)) = -1;
        elseif Edge2Node(edge_id,1)==start_node && Edge2Node(edge_id,2)==end_node % shared edge with same direction (if exists)
            Element2Edge(e,i_map(i)) = edge_id;
            Element2SignOfEdge(e,i_map(i)) = 1;
        end
    end
end
Edge2Node = Edge2Node(1:Num_Edges,:);

Edge2Element = zeros(Num_Edges,2); % edge shared by which two element
EdgeStartNodeLocalID = zeros(Num_Edges,2);
EdgeEndNodeLocalID = zeros(Num_Edges,2);
% ghost element ID is zero by default (the edge is only possessed by one element)
for e=1:Num_Elements
    % loop through the three edges of the tri element
    for i=1:3 
        edge_id = Element2Edge(e,i_map(i)); % global edge id
        edge_sign = Element2SignOfEdge(e,i_map(i));
        if edge_sign == 1 % +1 for counter-clockwise
            Edge2Element(edge_id,1) = e;
            EdgeStartNodeLocalID(edge_id,1) = l_map(i);
            EdgeEndNodeLocalID(edge_id,1) = k_map(i);
        elseif edge_sign == -1 % -1 for clockwise
            Edge2Element(edge_id,2) = e;
            EdgeStartNodeLocalID(edge_id,2) = k_map(i); % the opposite
            EdgeEndNodeLocalID(edge_id,2) = l_map(i);
        end
    end
end


%% Identify the source and bounaries
NodesOnPEC = findNodes(mesh,'region','Edge',[5, 6, 7, 8]);
NodesOnBoundary = findNodes(mesh,'region','Edge',[1, 2, 3, 4]);
EdgesIsOnPEC = zeros(Num_Edges,1); % flag for nodes with unknown coefficients
EdgesIsOnABC = zeros(Num_Edges,1); % flag for nodes with unknown coefficients
for edge_id = 1:Num_Edges
    if ismember(Edge2Node(edge_id,1),NodesOnPEC) && ismember(Edge2Node(edge_id,2),NodesOnPEC)
        EdgesIsOnPEC(edge_id)=1;
    end
    if ismember(Edge2Node(edge_id,1),NodesOnBoundary) && ismember(Edge2Node(edge_id,2),NodesOnBoundary)
        EdgesIsOnABC(edge_id)=1;
    end
end
NodeNearSource = findNodes(mesh,'nearest',[-1;-1]);
ElementofSource = findElements(mesh,'attached',NodeNearSource);
ElementofSource = ElementofSource(1);
%% Define coefficient matrices, vectors and fluxes
% basis only (Mass): [T]
Tez = zeros(Num_Elements,3,3); 
Thx = zeros(Num_Elements,3,3); 
Thy = zeros(Num_Elements,3,3); 
% ABC and PEC: [R]
Rez = zeros(Num_Elements,3,3); 
Rhx = zeros(Num_Elements,3,3); 
Rhy = zeros(Num_Elements,3,3); 
% basis & derivative of basis (Stiffness): [S]
Sx = zeros(Num_Elements,3,3); 
Sy = zeros(Num_Elements,3,3); 


%% Assembling the spatial matrix system.
xe = zeros(1,3);
ye = zeros(1,3);
ae = zeros(1,3);
be = zeros(1,3);
ce = zeros(1,3);
% element by element
for e=1:Num_Elements
    ne(1:3) = mesh.Elements(1:3,e); % the global node idx of node 1 in element e
    for i=1:3 % the global coordinate of node 1 in element e
        where_am_I = mesh.Nodes(:,ne(i)); 
        xe(i) = where_am_I(1);
        ye(i) = where_am_I(2);
    end

    eps_e = eps0*eps(e);
    mu_e = mu0*mu(e);
    sigma_e = 1*sigma(e);

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
            % basis only (Mass): [T]
            Tez(e,i,j) = eps_e*Delta_e/12*(1+(i==j)); 
            Thx(e,i,j) = mu_e*Delta_e/12*(1+(i==j));
            Thy(e,i,j) = mu_e*Delta_e/12*(1+(i==j)); 
            % ABC and PEC: [R]
            Rez(e,i,j) = sigma_e*Delta_e/12*(1+(i==j));
            Rhx(e,i,j) = sigma_e*Delta_e/12*(1+(i==j));
            Rhy(e,i,j) = sigma_e*Delta_e/12*(1+(i==j));
            % basis & derivative of basis (Stiffness): [S]
            Sx(e,i,j) = be(j)/6; 
            Sy(e,i,j) = ce(j)/6; 
        end
    end
end

%% Time-stepping and Solving Local Linear System
% ghost elements: the ones (non-existance) that provides solution values
% with '+' for the boundary elements
% ghost element boundary: the relation shape between electric field itself and magnetic field
% itself at nodes of boundary elements and ghost elements outside boundary
% Silver-Muller ABC boundary: the relationship between electric field at
% nodes of bounadry elements and magnetic field at nodes of ghost elements outside
% boundary

% step by step
for n=1:NumTimeSteps
    %% Initialization 
    t = n*dt;
    % flux and source terms should be renewed every time --> tricky !!!
    % interface/flux: {f}
    fez = zeros(Num_Elements,3,1); 
    fhx = zeros(Num_Elements,3,1); 
    fhy = zeros(Num_Elements,3,1); 
    % source: {l}
    lez = zeros(Num_Elements,3,1); 
    lhx = zeros(Num_Elements,3,1); 
    lhy = zeros(Num_Elements,3,1); 
    %% E update
    for e=1:Num_Elements
        %% Element information
        ne(1:3) = mesh.Elements(1:3,e); % the global node idx of node 1 in element e
        for i=1:3 % the global coordinate of node 1 in element e
            where_am_I = mesh.Nodes(:,ne(i)); 
            xe(i) = where_am_I(1);
            ye(i) = where_am_I(2);
        end

        eta_e = eta0*sqrt(mu(e)/eps(e));

        ae(1) = xe(2)*ye(3)-xe(3)*ye(2); ae(2) = xe(3)*ye(1)-xe(1)*ye(3); ae(3) = xe(1)*ye(2)-xe(2)*ye(1);
        be(1) = ye(2) - ye(3); be(2) = ye(3) - ye(1); be(3) = ye(1) - ye(2);
        ce(1) = xe(3) - xe(2); ce(2) = xe(1) - xe(3); ce(3) = xe(2) - xe(1);
        Delta_e = (be(1)*ce(2)-be(2)*ce(1))/2;
        %% Electric current source 
        if e == ElementofSource
            x_mid = sum(xe)/3;
            y_mid = sum(ye)/3;
            for i=1:3
                lez(e,i) = Delta_e*Jzimp(t)*Nel2D(x_mid,y_mid,ae(i),be(i),ce(i),Delta_e);      
            end
        end
        %% Coupling 
        for i=1:3 % loop through the three edges of the tri element
            edge_id = Element2Edge(e,i_map(i));
            edge_sign = Element2SignOfEdge(e,i_map(i));
            vecEdge = [xe(k_map(i))-xe(l_map(i)),ye(k_map(i))-ye(l_map(i))];
            edgeLength = norm(vecEdge);
            edgeDir = vecEdge/edgeLength;
            bn = edgeDir*Rotate90Mat;
            nx = bn(1);
            ny = bn(2);

            if  edge_sign == 1 % edge (as interface) in this element is positive edge
                this_local_node_id = [l_map(i),k_map(i)];
                neighborElement = Edge2Element(edge_id,2); % edge in neighboring element across this interface is negative edge
                neighbor_local_node_id = [EdgeStartNodeLocalID(edge_id,2),EdgeEndNodeLocalID(edge_id,2)];
            elseif edge_sign == -1 % on the contrary
                neighborElement = Edge2Element(edge_id,1);
                this_local_node_id = [k_map(i),l_map(i)];
                neighbor_local_node_id = [EdgeStartNodeLocalID(edge_id,1),EdgeEndNodeLocalID(edge_id,1)];
            end
           
            if neighborElement == 0 % ghost element
                assert(EdgesIsOnABC(edge_id) || EdgesIsOnPEC(edge_id),'Error in finding neighbor: this edge should be an interface')
                if EdgesIsOnPEC(edge_id) % PEC
                    for ii=1:2
                        ln_id_i = this_local_node_id(ii);
                        ffhx = 0; ffhy = 0;
                        for jj=1:2
                            ln_id_j = this_local_node_id(jj);
                            ff = edgeLength/6*(1+(ln_id_i==ln_id_j));
                            hx_coeff_plus = hx_coeff1(e,ln_id_j);
                            hy_coeff_plus = hy_coeff1(e,ln_id_j);
                            ffhx = ffhx + ff*(hx_coeff_plus - hx_coeff1(e,ln_id_j));
                            ffhy = ffhy + ff*(hy_coeff_plus - hy_coeff1(e,ln_id_j));
                        end
                        fez(e,ln_id_i) = fez(e,ln_id_i) - (ny*ffhx-nx*ffhy)/2; 
                    end
                end
                if EdgesIsOnABC(edge_id) % ABC
                    for ii=1:2
                        ln_id_i = this_local_node_id(ii);
                        ffhx = 0; ffhy = 0;
                        for jj=1:2
                            ln_id_j = this_local_node_id(jj);
                            ff = edgeLength/6*(1+(ln_id_i==ln_id_j));
                            hx_coeff_plus = 1/eta_e*ny*ez_coeff1(e,ln_id_j);
                            hy_coeff_plus = -1/eta_e*nx*ez_coeff1(e,ln_id_j);
                            ffhx = ffhx + ff*(hx_coeff_plus - hx_coeff1(e,ln_id_j));
                            ffhy = ffhy + ff*(hy_coeff_plus - hy_coeff1(e,ln_id_j));
                        end
                        fez(e,ln_id_i) = fez(e,ln_id_i) - (ny*ffhx-nx*ffhy)/2; 
                    end
                end
            else % is the element on the other side of the interface 
                % element interface
                for ii=1:2 
                    ln_id_i = this_local_node_id(ii); % q'
                    ffhx = 0; ffhy = 0;
                    for jj=1:2 
                        neighbor_ln_id_j = neighbor_local_node_id(jj); % p
                        ln_id_j = this_local_node_id(jj); % q
                        assert(mesh.Elements(neighbor_ln_id_j,neighborElement)==mesh.Elements(ln_id_j,e), 'Error in flux coupling: not the same node');

                        ff = edgeLength/6*(1+(ln_id_i==ln_id_j));
                        hx_coeff_plus = hx_coeff1(neighborElement,neighbor_ln_id_j);
                        hy_coeff_plus = hy_coeff1(neighborElement,neighbor_ln_id_j);
                        ffhx = ffhx + ff*(hx_coeff_plus - hx_coeff1(e,ln_id_j));
                        ffhy = ffhy + ff*(hy_coeff_plus - hy_coeff1(e,ln_id_j));
                    end
                    fez(e,ln_id_i) = fez(e,ln_id_i) - (ny*ffhx-nx*ffhy)/2; 
                end
            end
        end
    end
    
    for e=1:Num_Elements
        Teze = squeeze(Tez(e,:,:)); 
        Reze = squeeze(Rez(e,:,:)); 
        Sxe = squeeze(Sx(e,:,:)); 
        Sye = squeeze(Sy(e,:,:)); 
        feze = squeeze(fez(e,:))'; 
        leze = squeeze(lez(e,:))'; 
        
        Keze = 1/dt*Teze+1/2*Reze;
        beze = (1/dt*Teze-1/2*Reze)*ez_coeff1(e,:)'+Sxe*hy_coeff1(e,:)'-Sye*hx_coeff1(e,:)'+feze+leze;
        ez_coeff2(e,:) = Keze\beze;
    end
    
    
    %% H update
    for e=1:Num_Elements
        
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

        ae(1) = xe(2)*ye(3)-xe(3)*ye(2); ae(2) = xe(3)*ye(1)-xe(1)*ye(3); ae(3) = xe(1)*ye(2)-xe(2)*ye(1);
        be(1) = ye(2) - ye(3); be(2) = ye(3) - ye(1); be(3) = ye(1) - ye(2);
        ce(1) = xe(3) - xe(2); ce(2) = xe(1) - xe(3); ce(3) = xe(2) - xe(1);
        Delta_e = (be(1)*ce(2)-be(2)*ce(1))/2;
       
        %% Coupling
        for i=1:3 % loop through the three edges of the tri element
            edge_id = Element2Edge(e,i_map(i));
            edge_sign = Element2SignOfEdge(e,i_map(i));
            vecEdge = [xe(k_map(i))-xe(l_map(i)),ye(k_map(i))-ye(l_map(i))];
            edgeLength = norm(vecEdge);
            edgeDir = vecEdge/edgeLength;
            bn = edgeDir*Rotate90Mat;
            nx = bn(1);
            ny = bn(2);

            if  edge_sign == 1 % edge (as interface) in this element is positive edge
                this_local_node_id = [l_map(i),k_map(i)]; % [EdgeStartNodeLocalID(edge_id,1),EdgeEndNodeLocalID(edge_id,1)]
                neighborElement = Edge2Element(edge_id,2); % edge in neighboring element across this interface is negative edge
                neighbor_local_node_id = [EdgeStartNodeLocalID(edge_id,2),EdgeEndNodeLocalID(edge_id,2)];
            elseif edge_sign == -1 % on the contrary
                neighborElement = Edge2Element(edge_id,1);
                this_local_node_id = [k_map(i),l_map(i)]; % [EdgeStartNodeLocalID(edge_id,2),EdgeEndNodeLocalID(edge_id,2)]
                neighbor_local_node_id = [EdgeStartNodeLocalID(edge_id,1),EdgeEndNodeLocalID(edge_id,1)];
            end

            if neighborElement == 0 % ghost element
                assert(EdgesIsOnABC(edge_id) || EdgesIsOnPEC(edge_id),'Error in finding neighbor: this edge should be an interface')
                if EdgesIsOnPEC(edge_id) % PEC
                    for ii=1:2
                        ln_id_i = this_local_node_id(ii);
                        ffez = 0;
                        for jj=1:2
                            ln_id_j = this_local_node_id(jj);
                            ff = edgeLength/6*(1+(ln_id_i==ln_id_j));
                            ez_coeff_plus = - ez_coeff2(e,ln_id_j);
                            ffez = ffez + ff*(ez_coeff_plus - ez_coeff2(e,ln_id_j));
                        end
                        fhx(e,ln_id_i) = fhx(e,ln_id_i) - ny*ffez/2; 
                        fhy(e,ln_id_i) = fhy(e,ln_id_i) + nx*ffez/2; 
                    end
                end
                if EdgesIsOnABC(edge_id) % ABC
                    for ii=1:2
                        ln_id_i = this_local_node_id(ii);
                        ffez = 0;
                        for jj=1:2
                            ln_id_j = this_local_node_id(jj);
                            ff = edgeLength/6*(1+(ln_id_i==ln_id_j));
                            ez_coeff_plus = eta_e*(ny*hx_coeff1(e,ln_id_j)-nx*hy_coeff1(e,ln_id_j));
                            ffez = ffez + ff*(ez_coeff_plus - ez_coeff2(e,ln_id_j));
                        end
                        fhx(e,ln_id_i) = fhx(e,ln_id_i) - ny*ffez/2; 
                        fhy(e,ln_id_i) = fhy(e,ln_id_i) + nx*ffez/2; 
                    end
                end
            else % is the element on the other side of the interface 
                % element interface
                for ii=1:2 
                    ln_id_i = this_local_node_id(ii); % q'
                    ffez = 0;
                    for jj=1:2 
                        neighbor_ln_id_j = neighbor_local_node_id(jj); % p
                        ln_id_j = this_local_node_id(jj); % q
                        assert(mesh.Elements(neighbor_ln_id_j,neighborElement)==mesh.Elements(ln_id_j,e), 'Error in flux coupling: not the same node');

                        ff = edgeLength/6*(1+(ln_id_i==ln_id_j));
                        ez_coeff_plus = ez_coeff2(neighborElement,neighbor_ln_id_j);
                        ffez = ffez + ff*(ez_coeff_plus - ez_coeff2(e,ln_id_j));
                    end
                    fhx(e,ln_id_i) = fhx(e,ln_id_i) - ny*ffez/2; 
                    fhy(e,ln_id_i) = fhy(e,ln_id_i) + nx*ffez/2; 
                end
            end
        end
    end
    
    for e=1:Num_Elements
        Thxe = squeeze(Thx(e,:,:)); Thye = squeeze(Thy(e,:,:));
        Rhxe = squeeze(Rhx(e,:,:));  Rhye = squeeze(Rhy(e,:,:));
        Sxe = squeeze(Sx(e,:,:));  Sye = squeeze(Sy(e,:,:)); 
        fhxe = squeeze(fhx(e,:))';  fhye = squeeze(fhy(e,:))'; 
        lhxe = squeeze(lhx(e,:))';  lhye = squeeze(lhy(e,:))'; 
        
        Khxe = 1/dt*Thxe+1/2*Rhxe; 
        bhxe = (1/dt*Thxe-1/2*Rhxe)*hx_coeff1(e,:)'-Sye*ez_coeff2(e,:)'+fhxe+lhxe;
        hx_coeff2(e,:) = Khxe\bhxe;
        
        Khye = 1/dt*Thye+1/2*Rhye;
        bhye = (1/dt*Thye-1/2*Rhye)*hy_coeff1(e,:)'+Sxe*ez_coeff2(e,:)'+fhye+lhye;
        hy_coeff2(e,:) = Khye\bhye;
    end
    
    %% Initiate next iteration
    ez_coeff1 = ez_coeff2;
    hx_coeff1 = hx_coeff2;
    hy_coeff1 = hy_coeff2;
    
    %% Visualization
    for e=1:Num_Elements
        ne(1:3) = mesh.Elements(1:3,e); % the global node idx of node 1 in element e
        for i=1:3
            ez_coeff_global(ne(i)) = ez_coeff1(e,i);
        end
    end
    if use_pdeplot
        figure(5);
        pdeplot(model, 'XYData', abs(real(ez_coeff_global)));
        colormap summer;colorbar;axis image;xlabel('Y (m)');ylabel('X (m)');set(gca,'fontsize',18);
        title(['$E_z^{tot}$', ' at ', 'time = ',num2str(round(n*dt,9)),' s',],'interpreter','latex');
        drawnow; pause(0.2);
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

                ae(1) = xe(2)*ye(3)-xe(3)*ye(2); ae(2) = xe(3)*ye(1)-xe(1)*ye(3); ae(3) = xe(1)*ye(2)-xe(2)*ye(1);
                be(1) = ye(2) - ye(3); be(2) = ye(3) - ye(1); be(3) = ye(1) - ye(2);
                ce(1) = xe(3) - xe(2); ce(2) = xe(1) - xe(3); ce(3) = xe(2) - xe(1);
                Delta_e = (be(1)*ce(2)-be(2)*ce(1))/2;
                for q=1:N
                    for p=1:N
                        if ~GridInterpolated(q,p)
                            [in,on]=inpolygon(X(p),Y(q),xe,ye);
                            if (in || on) %indicating if inside the element or on the edge of the element
                                for i=1:3
                                    EzGrid(q,p)=EzGrid(q,p)+Nel2D(X(p),Y(q),ae(i),be(i),ce(i),Delta_e)*ez_coeff1(e,i);
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