% 2D EDGE-BASED (VECTOR) HYBRID FETD-FDTD METHOD, TM FIELD
clear all; close all; clc; tic;
%% settings
disp('start setting basic parameters ... ')
N = 150; % interpolation size
Hmax=0.12; % maximum element size

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
NumTimeSteps = 1000;
plot_every_nsteps = 20;
% FETD region: (a+d1) x (a+d1)
a = lambda*0.5; % radius of circle object
d1 = lambda*0.5; % distance between obj and truncation box
d2 = lambda*2; % FDTD thickness
d3 = lambda; % PML thickness

% choose what to plot and how to plot
plot_logistics = false; % if to plot FEM matrices and meshes
save_fig = true; % if to save fig
use_fdtd_update_with_pml = true; % if to use FDTD with well-posed PML
apply_lpf = true; % if to solve late time instability according to (Hwang 1999)

% single freq electric current Jx
J_dir = [1,0];
J_imp_amp = 1; tau_p = 3/omega;

%% macros and lambda expressions
disp('start defining macros and lambda expressions ... ')
% triangular basis functions
Nel2D = @(x,y,ael,bel,cel,Delta_e) 1/(2*Delta_e)*(ael+bel*x+cel*y); 
NablaNel2D = @(ael,bel,cel,Delta_e) 1/(2*Delta_e).*[bel,cel];
Nelk2D = @(x,y,ael,bel,cel,aek,bek,cek,le_lk,Delta_e) (Nel2D(x,y,ael,bel,cel,Delta_e)*NablaNel2D(aek,bek,cek,Delta_e)-Nel2D(x,y,aek,bek,cek,Delta_e)*NablaNel2D(ael,bel,cel,Delta_e))*le_lk;
CurlOfNelk2D = @(bel,cel,bek,cek,le_lk,Delta_e) le_lk/(2*Delta_e^2)*(bel*cek-bek*cel);

% rectangular basis functions
N1 = @(x,y,xedge,dx) (1-(x-xedge)/dx).*[0,1];
N2 = @(x,y,yedge,dy) (1-(yedge-y)/dy).*[1,0];
N3 = @(x,y,xedge,dx) (1-(xedge-x)/dx).*[0,1];
N4 = @(x,y,yedge,dy) (1-(y-yedge)/dy).*[1,0];

GaussWeightsTri3 = [1/3,1/3,1/3];
GaussQuadTri3 = @(GaussSamples) GaussSamples*GaussWeightsTri3';

Rotate90Mat = [cosd(90) -sind(90); sind(90) cosd(90)];

% tri edge defined counter clock wise
tri_i_map = [3,1,2]; % from loop index to triangular local edge id
tri_l_map = [1,2,3]; % from loop index to triangular local start node id
tri_k_map = [2,3,1]; % from loop index to triangular local end node id

% rect edge defined clock wise
rect_i_map = [1,2,3,4]; % from loop index to rectangular local edge id
rect_l_map = [1,4,3,2]; % from loop index to rectangular local start node id
rect_k_map = [4,3,2,1]; % from loop index to rectangular local end node id
rect_i_sign = [1,1,-1,-1]; % from loop index to rectangular local edge sign

Jimp = @(t) J_imp_amp*(1-exp(-(t)/tau_p)).*sin(omega*t);
dJimp_dt = @(t) omega*J_imp_amp*(1-exp(-(t)/tau_p)).*cos(omega*t)+J_imp_amp*(exp(-(t)/tau_p)/tau_p).*sin(omega*t);

%% Geometries in FETD Region
disp('start generating geometries in FETD region ... ')
model = createpde;
Circle = [1 0 0 a]';
Rect = [3 4 -a-d1 -a-d1 a+d1 a+d1 -a-d1 a+d1 a+d1 -a-d1]';
Circle = [Circle;zeros(length(Rect) - length(Circle),1)];
gd = [Circle, Rect];
ns = char('Circle','Rect');
ns = ns';
sf = 'Rect-Circle';
g= decsg(gd,sf,ns);
geometryFromEdges(model,g);
if plot_logistics
    figure; hold on;
    pdegplot(model,'FaceLabels','on','EdgeLabels','on')
    xlabel('x');
    xlabel('y');
    hold off;
end
%% Triangular Mesh in FETD Region
disp('start generating meshes in FETD region ... ')
% Hmax = 2*(d3+d2+d1+a)/ceil(2*(d3+d2+d1+a)/Hmax); % refine Hmax
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',Hmax);
NodesOnUpInterface = findNodes(mesh,'region','Edge',1);
NodesOnBottomInterface = findNodes(mesh,'region','Edge',4);
NodesOnLeftInterface = findNodes(mesh,'region','Edge',3);
NodesOnRightInterface = findNodes(mesh,'region','Edge',2);

assert(length(NodesOnBottomInterface)==length(NodesOnUpInterface),'Error: horizontal discretization do not fulfill our simple assumption');
assert(length(NodesOnLeftInterface)==length(NodesOnRightInterface),'Error: vertical discretization do not fulfill our simple assumption');

Num_Tri_Element_X_Interface = length(NodesOnBottomInterface)-1;
Num_Tri_Element_Y_Interface = length(NodesOnLeftInterface)-1;

Mesh_Nodes_Tri = mesh.Nodes;
Mesh_Elements_Tri = mesh.Elements;
Num_Nodes_Tri = size(Mesh_Nodes_Tri,2);
Num_Elements_Tri = size(Mesh_Elements_Tri,2);

TriNodesOnPEC = findNodes(mesh,'region','Edge',[5, 6, 7, 8]);
TriNodesOnBoundary = findNodes(mesh,'region','Edge',[1, 2, 3, 4]);
TriNodesNearSource = findNodes(mesh,'nearest',[-a;-a]);
TriElementofSource = findElements(mesh,'attached',TriNodesNearSource);
TriElementofSource = TriElementofSource(1);

TriElementIsAdjToBoundary = zeros(1,Num_Elements_Tri);
dx_List = zeros(1,Num_Tri_Element_X_Interface); 
dy_List = zeros(1,Num_Tri_Element_Y_Interface);
num_dx = 0;
num_dy = 0;
for tri_e = 1:Num_Elements_Tri
    ne_tri(1:3) = Mesh_Elements_Tri(1:3,tri_e);
    for i=1:3  
        if ismember(ne_tri(tri_l_map(i)),TriNodesOnBoundary) && ismember(ne_tri(tri_k_map(i)),TriNodesOnBoundary)
            TriElementIsAdjToBoundary(tri_e)=1;
        end
        
        % check if the edges of tri elements are equal-length at the interface 
        if ismember(ne_tri(tri_l_map(i)),NodesOnBottomInterface) && ismember(ne_tri(tri_k_map(i)),NodesOnBottomInterface)
            num_dx = num_dx + 1;
            dx_List(num_dx)=norm(Mesh_Nodes_Tri(:,ne_tri(tri_k_map(i)))-Mesh_Nodes_Tri(:,ne_tri(tri_l_map(i))));
        end
        if ismember(ne_tri(tri_l_map(i)),NodesOnLeftInterface) && ismember(ne_tri(tri_k_map(i)),NodesOnLeftInterface)
            num_dy = num_dy + 1;
            dy_List(num_dy)=norm(Mesh_Nodes_Tri(:,ne_tri(tri_k_map(i)))-Mesh_Nodes_Tri(:,ne_tri(tri_l_map(i))));
        end
    end
end
%% Rectangular Mesh in FDTD Region
disp('start generating meshes in FDTD region ... ')
% must generate conformal mesh at the interface
DX = 2*(a+d1)/Num_Tri_Element_X_Interface;
DY = 2*(a+d1)/Num_Tri_Element_Y_Interface;
dt = 1/c/sqrt(1/DX^2+1/DY^2)/2; % choose time step size according to page 395 of text book
assert(DX==DY,'Error: dx~=dy');
assert(norm(dx_List-DX)<1e-12 && norm(dy_List-DY)<1e-12,'Error: mesh is not conformal')

% refine PML thickness d3 and FDTD region size d2
d3 = floor(d3/DX)*DX;
d2 = floor(d2/DX)*DX;
xx = [-(d3+d2+d1+a):DX:-(d1+a+DX),-(d1+a):DX:(d1+a),(d1+a+DX):DX:(d3+d2+d1+a)];
yy = [-(d3+d2+d1+a):DY:-(d1+a+DY),-(d1+a):DY:(d1+a),(d1+a+DY):DY:(d3+d2+d1+a)];
[Mesh_Nodes_Rect,Mesh_Elements_Rect,RectNodeIsOnBoundary,RectElementIsInPML,RectElementIsAdjToBoundary,RectElement2Grid2D] = RectangularMeshForFDTDRegion(xx,yy,d1+a,a+d1+d2);
Num_Nodes_Rect = size(Mesh_Nodes_Rect,2);
Num_Elements_Rect = size(Mesh_Elements_Rect,2);

%% Edge numbering and its connection to nodes and elements. 
disp('start generating edge numbering and connection to nodes and elements ... ')
% Some philosophies about DDM:
% 1. Make the numbering of edges with unknowns global
% 2. Keep the element numbering and node numbering local

Num_Edges = 0; % total number of edges && edge id
TriElement2Edge = zeros(Num_Elements_Tri,3);
RectElement2Edge = zeros(Num_Elements_Rect,4);

Edge2TriNode = zeros(3*Num_Elements_Tri+4*Num_Elements_Rect,2);
Edge2RectNode = zeros(3*Num_Elements_Tri+4*Num_Elements_Rect,2);

TriElement2SignOfEdge = zeros(Num_Elements_Tri,3);
RectElement2SignOfEdge = zeros(Num_Elements_Rect,4);

for e=1:Num_Elements_Tri
    ne_tri(1:3) = Mesh_Elements_Tri(1:3,e);
    for i=1:3
        start_node = ne_tri(tri_l_map(i));
        end_node = ne_tri(tri_k_map(i));
        edge_id = find( ( Edge2TriNode(:,1)==start_node & Edge2TriNode(:,2)==end_node ) | ( Edge2TriNode(:,2)==start_node & Edge2TriNode(:,1)==end_node ), 1);
        
        if isempty(edge_id) % new edge
            % first add it to edge list
            Num_Edges = Num_Edges + 1;
            Edge2TriNode(Num_Edges,1) = start_node;
            Edge2TriNode(Num_Edges,2) = end_node;
            % second add the edge to the element-edge list
            TriElement2Edge(e,tri_i_map(i)) = Num_Edges;
            TriElement2SignOfEdge(e,tri_i_map(i)) = 1;
        elseif Edge2TriNode(edge_id,2)==start_node && Edge2TriNode(edge_id,1)==end_node % shared edge with opposite direction
            TriElement2Edge(e,tri_i_map(i)) = edge_id;
            TriElement2SignOfEdge(e,tri_i_map(i)) = -1;
        elseif Edge2TriNode(edge_id,1)==start_node && Edge2TriNode(edge_id,2)==end_node % shared edge with same direction (if exists)
            TriElement2Edge(e,tri_i_map(i)) = edge_id;
            TriElement2SignOfEdge(e,tri_i_map(i)) = 1;
        end
    end
end
Num_Edges_Tri = Num_Edges;
for e=1:Num_Elements_Rect
    ne_rect(1:4) = Mesh_Elements_Rect(1:4,e);
    for i=1:4
        start_node = ne_rect(rect_l_map(i));
        end_node = ne_rect(rect_k_map(i));
        if ~(RectNodeIsOnBoundary(start_node) && RectNodeIsOnBoundary(end_node)) % this is an edge belong to rectangular element region
            edge_id = find( ( Edge2RectNode(:,1)==start_node & Edge2RectNode(:,2)==end_node ) | ( Edge2RectNode(:,2)==start_node & Edge2RectNode(:,1)==end_node ), 1);
            
            if isempty(edge_id) % new edge
                % first add it to edge list
                Num_Edges = Num_Edges + 1;
                Edge2RectNode(Num_Edges,1) = start_node;
                Edge2RectNode(Num_Edges,2) = end_node;
                % second add the edge to the element-edge list
                RectElement2Edge(e,rect_i_map(i)) = Num_Edges;
                RectElement2SignOfEdge(e,rect_i_map(i)) = rect_i_sign(i);
            elseif Edge2RectNode(edge_id,2)==start_node && Edge2RectNode(edge_id,1)==end_node % shared edge with opposite direction
                RectElement2Edge(e,rect_i_map(i)) = edge_id;
                RectElement2SignOfEdge(e,rect_i_map(i)) = rect_i_sign(i);
            elseif Edge2RectNode(edge_id,1)==start_node && Edge2RectNode(edge_id,2)==end_node % shared edge with same direction (if exists)
                RectElement2Edge(e,rect_i_map(i)) = edge_id;
                RectElement2SignOfEdge(e,rect_i_map(i)) = rect_i_sign(i);
            end
        end
    end
end
Num_Edges_Rect = Num_Edges - Num_Edges_Tri;
Edge2TriNode = Edge2TriNode(1:Num_Edges,:);
Edge2RectNode = Edge2RectNode(1:Num_Edges,:);

% Flags for special edges and elements
EdgesOnPEC = [];
EdgesIsOnBoundary = zeros(Num_Edges,1);
for tri_e = 1:Num_Elements_Tri
    for edge_id=TriElement2Edge(tri_e,:)
        if ismember(Edge2TriNode(edge_id,1),TriNodesOnPEC) && ismember(Edge2TriNode(edge_id,2),TriNodesOnPEC)
            EdgesOnPEC = [EdgesOnPEC,edge_id];
        end
        if ismember(Edge2TriNode(edge_id,1),TriNodesOnBoundary) && ismember(Edge2TriNode(edge_id,2),TriNodesOnBoundary)
            EdgesIsOnBoundary(edge_id)=1;
        end
    end
end

% Find edges of rectangular elements that are respectively associated with
% Ex and Ey
Ex_List = [];
Ey_List = [];
for e=1:Num_Elements_Rect
    if ~ismember(RectElement2Edge(e,4),Ex_List)
        Ex_List = [Ex_List,RectElement2Edge(e,4)];
    end
    if ~ismember(RectElement2Edge(e,2),Ex_List)
        Ex_List = [Ex_List,RectElement2Edge(e,2)];
    end
    if ~ismember(RectElement2Edge(e,1),Ey_List)
        Ey_List = [Ey_List,RectElement2Edge(e,1)];
    end
    if ~ismember(RectElement2Edge(e,3),Ey_List)
        Ey_List = [Ey_List,RectElement2Edge(e,3)];
    end
end

%% Elements that are across the interface
disp('start identifying elements across the interface ... ')
% Identify the rectangular elements in 'Mesh_Elements_Rect' conformally
% connected to the edge of the triangular elements in 'Mesh_Elements_Tri'
% The edge is uniquely defined and where the basis function and coefficient
% locates.
RectElement2AdjTriEdge = zeros(Num_Elements_Rect,1);
RectElement2AdjTriElement = zeros(Num_Elements_Rect,1); 
for edge_id = 1:Num_Edges_Tri  % loop through edges
    if EdgesIsOnBoundary(edge_id)
        for tri_e = 1:Num_Elements_Tri % loop through tri elements
            if TriElementIsAdjToBoundary(tri_e)
                if ismember(edge_id,TriElement2Edge(tri_e,:)) % tri_e contains this interface edge
                    for rect_e = 1:Num_Elements_Rect % loop through rect elements
                        if RectElementIsAdjToBoundary(rect_e)
                            rect_local_edge_id = find(RectElement2Edge(rect_e,:)==0,1);
                            rect_local_edge_idx = find(rect_i_map==rect_local_edge_id,1);
                            ne_rect(1:4) = Mesh_Elements_Rect(1:4,rect_e);
                            edge_tri = [Mesh_Nodes_Tri(:,Edge2TriNode(edge_id,2)),Mesh_Nodes_Tri(:,Edge2TriNode(edge_id,1))];
                            edge_rect = [Mesh_Nodes_Rect(:,ne_rect(rect_k_map(rect_local_edge_idx))),Mesh_Nodes_Rect(:,ne_rect(rect_l_map(rect_local_edge_idx)))];

                            if isApproximatelyEqual(edge_rect,edge_tri,1e-12) % approximately equal
                                disp(['Rect element ',num2str(rect_e),' is connected to tri element ', num2str(tri_e),' through edge ',num2str(edge_id)]);
                                RectElement2AdjTriEdge(rect_e) = edge_id;
                                RectElement2AdjTriElement(rect_e) = tri_e;
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Visualization of hybrid grids
disp('start generating figure of hybrid grids ... ')
if plot_logistics
    figure; hold on;
    pdemesh(mesh,'ElementLabels','off')
    for e=1:Num_Elements_Rect
        ne_rect(1:4) = Mesh_Elements_Rect(1:4,e);
        where_am_I = Mesh_Nodes_Rect(:,ne_rect); 
        xe_rect = where_am_I(1,:);
        ye_rect = where_am_I(2,:);
        if RectElementIsInPML(e)~=0
            plot([xe_rect xe_rect(1)], [ye_rect ye_rect(1)], '-g');
            text(mean(xe_rect),mean(ye_rect),num2str(RectElementIsInPML(e)));
        else
            plot([xe_rect xe_rect(1)], [ye_rect ye_rect(1)], '-k');
        end
        
    end
    for e=1:Num_Elements_Tri
        if TriElementIsAdjToBoundary(e)
            for i=1:3
                edge_id = TriElement2Edge(e,tri_i_map(i));
                plot(Mesh_Nodes_Tri(1,Edge2TriNode(edge_id,:)),Mesh_Nodes_Tri(2,Edge2TriNode(edge_id,:)),'-b','linewidth',1.5);
            end
        end
    end
    for e=1:Num_Elements_Rect
        if RectElementIsAdjToBoundary(e)
            for i=1:4
                edge_id = RectElement2Edge(e,rect_i_map(i));
                if edge_id~=0
                    plot(Mesh_Nodes_Rect(1,Edge2RectNode(edge_id,:)),Mesh_Nodes_Rect(2,Edge2RectNode(edge_id,:)),'-k','linewidth',1.5);
                end
            end
        end
    end
    for edge_id = 1:Num_Edges
        if EdgesIsOnBoundary(edge_id) 
            plot(Mesh_Nodes_Tri(1,Edge2TriNode(edge_id,:)),Mesh_Nodes_Tri(2,Edge2TriNode(edge_id,:)),'-b','linewidth',2.5);
        end
    end
    xlabel('x');
    xlabel('y');
    set(gca,'fontsize',24);
	drawnow;
	if save_fig
		set(gcf,'Position',[0,0,1024,1024]); 
		saveas(gca, ['Hybrid-FETD-FDTD Mesh and Domain','.png']);
	end
    hold off;
end


%% Materials and PML settings
disp('start defining materials and generating PML for FDTD region ... ')
eps_tri = ones(Num_Elements_Tri,1); % relative permittivity of the tri elements
mu_tri = ones(Num_Elements_Tri,1); % relative permeability of the tri elements
eps_rect = ones(Num_Elements_Rect,1); % relative permittivity of the rect elements
mu_rect = ones(Num_Elements_Rect,1); % relative permeability of the rect elements
sigma_x = zeros(Num_Edges,1);
sigma_y = zeros(Num_Edges,1);
beta_x = mu0/dt*ones(Num_Edges,1);
beta_y = mu0/dt*ones(Num_Edges,1);
alpha_x = mu0/dt*ones(Num_Edges,1);
alpha_y = mu0/dt*ones(Num_Edges,1);

if use_fdtd_update_with_pml
    % PML of length NPML at outer grids
    NPML = sqrt(length(find(a==1)))+1; % the number of edges for one direction inside the PML region
    mPML = 2; % m = 2 or 3 is a good choice
    R0 = 1e-12; % magnitude of reflection coefficient
    sigma_max = -(mPML+1)/(2*NPML*DX/eta0)*log(R0); % equation on page of 416 textbook
    d2d1a=d2+d1+a;
    for rect_e=1:Num_Elements_Rect
        if RectElementIsInPML(rect_e)~=0
            for i=1:4
                edge_id = RectElement2Edge(rect_e,rect_i_map(i));
                x = (Mesh_Nodes_Rect(1,Edge2RectNode(edge_id,1))+Mesh_Nodes_Rect(1,Edge2RectNode(edge_id,2)))/2;
                y = (Mesh_Nodes_Rect(2,Edge2RectNode(edge_id,1))+Mesh_Nodes_Rect(2,Edge2RectNode(edge_id,2)))/2;
                
                beta_x(edge_id) = mu0*mu_rect(rect_e)/dt;
                beta_y(edge_id) = mu0*mu_rect(rect_e)/dt;
                alpha_x(edge_id) = mu0*mu_rect(rect_e)/dt;
                alpha_y(edge_id) = mu0*mu_rect(rect_e)/dt;
                
                if RectElementIsInPML(rect_e)==1
                    if ismember(edge_id,Ey_List)
                        sigma_x(edge_id) = sigma_max*((-d2d1a-x)/d3)^mPML;
                    elseif ismember(edge_id,Ex_List)
                        sigma_y(edge_id) = sigma_max*((-d2d1a+y)/d3)^mPML;
                    else
                        ME = MException('Error in edge identification: neither along x nor along y');
                        throw(ME);
                    end
                elseif RectElementIsInPML(rect_e)==2
                    if ismember(edge_id,Ey_List)
                        sigma_x(edge_id) = 0;
                    elseif ismember(edge_id,Ex_List)
                        sigma_y(edge_id) = sigma_max*((-d2d1a+y)/d3)^mPML;
                    else
                        ME = MException('Error in edge identification: neither along x nor along y');
                        throw(ME);
                    end
                elseif RectElementIsInPML(rect_e)==3
                    if ismember(edge_id,Ey_List)
                        sigma_x(edge_id) = sigma_max*((-d2d1a+x)/d3)^mPML;
                    elseif ismember(edge_id,Ex_List)
                        sigma_y(edge_id) = sigma_max*((-d2d1a+y)/d3)^mPML;
                    else
                        ME = MException('Error in edge identification: neither along x nor along y');
                        throw(ME);
                    end
                elseif RectElementIsInPML(rect_e)==4
                    if ismember(edge_id,Ey_List)
                        sigma_x(edge_id) = sigma_max*((-d2d1a+x)/d3)^mPML;
                    elseif ismember(edge_id,Ex_List)
                        sigma_y(edge_id) = 0;
                    else
                        ME = MException('Error in edge identification: neither along x nor along y');
                        throw(ME);
                    end
                elseif RectElementIsInPML(rect_e)==5
                    if ismember(edge_id,Ey_List)
                        sigma_x(edge_id) = sigma_max*((-d2d1a+x)/d3)^mPML;
                    elseif ismember(edge_id,Ex_List)
                        sigma_y(edge_id) = sigma_max*((-d2d1a-y)/d3)^mPML;
                    else
                        ME = MException('Error in edge identification: neither along x nor along y');
                        throw(ME);
                    end
                elseif RectElementIsInPML(rect_e)==6
                    if ismember(edge_id,Ey_List)
                        sigma_x(edge_id) = 0;
                    elseif ismember(edge_id,Ex_List)
                        sigma_y(edge_id) = sigma_max*((-d2d1a-y)/d3)^mPML;
                    else
                        ME = MException('Error in edge identification: neither along x nor along y');
                        throw(ME);
                    end
                elseif RectElementIsInPML(rect_e)==7
                    if ismember(edge_id,Ey_List)
                        sigma_x(edge_id) = sigma_max*((-d2d1a-x)/d3)^mPML;
                    elseif ismember(edge_id,Ex_List)
                        sigma_y(edge_id) = sigma_max*((-d2d1a-y)/d3)^mPML;
                    else
                        ME = MException('Error in edge identification: neither along x nor along y');
                        throw(ME);
                    end
                elseif RectElementIsInPML(rect_e)==8
                    if ismember(edge_id,Ey_List)
                        sigma_x(edge_id) = sigma_max*((-d2d1a-x)/d3)^mPML;
                    elseif ismember(edge_id,Ex_List)
                        sigma_y(edge_id) = 0;
                    else
                        ME = MException('Error in edge identification: neither along x nor along y');
                        throw(ME);
                    end
                end
            end
        end
    end

    beta_x = beta_x + sigma_x/2; 
    beta_y = beta_y + sigma_y/2;
    alpha_x = alpha_x - sigma_x/2;
    alpha_y = alpha_y - sigma_y/2;
end

%% FDTD Grid - FEM Element and Edge
disp('start preparing for FDTD leap frog fashion ... ')
% (1+2) matrices for Ey, (1+2) matrices for Ex, 2 matrices for Hz
% the '1' is the mapping matrix from FDTD grid to global edge
% the '2' are the value matrices for the adjacent two time steps
% all FDTD matrices are named in CAPITAL form
if use_fdtd_update_with_pml
    FDTDGrid_Size = length(xx);
    ExMap = zeros(FDTDGrid_Size-1,FDTDGrid_Size);
	ExSign = ones(FDTDGrid_Size-1,FDTDGrid_Size);
    Ex2 = zeros(FDTDGrid_Size-1,FDTDGrid_Size);
    Ex1 = zeros(FDTDGrid_Size-1,FDTDGrid_Size);

    EyMap = zeros(FDTDGrid_Size,FDTDGrid_Size-1);
	EySign = ones(FDTDGrid_Size,FDTDGrid_Size-1);
    Ey2 = zeros(FDTDGrid_Size,FDTDGrid_Size-1);
    Ey1 = zeros(FDTDGrid_Size,FDTDGrid_Size-1);

    Hzx2 = zeros(FDTDGrid_Size-1,FDTDGrid_Size-1);
    Hzx1 = zeros(FDTDGrid_Size-1,FDTDGrid_Size-1);
    Hzy2 = zeros(FDTDGrid_Size-1,FDTDGrid_Size-1);
    Hzy1 = zeros(FDTDGrid_Size-1,FDTDGrid_Size-1);

    BETAx = ones(FDTDGrid_Size,FDTDGrid_Size-1); 
    BETAy = ones(FDTDGrid_Size-1,FDTDGrid_Size);
    ALPHAx = ones(FDTDGrid_Size,FDTDGrid_Size-1);
    ALPHAy = ones(FDTDGrid_Size-1,FDTDGrid_Size);

    for rect_e=1:Num_Elements_Rect
        grid_edge_id = zeros(1,4);
		grid_edge_sign = zeros(1,4);
        ne_rect(1:4) = Mesh_Elements_Rect(1:4,rect_e);
        where_am_I = Mesh_Nodes_Rect(:,ne_rect); 
        xe_rect = where_am_I(1,:);
        ye_rect = where_am_I(2,:);
        dx = xe_rect(2)-xe_rect(1);
        dy = ye_rect(4)-ye_rect(1);
        eps_e = eps0*eps_rect(rect_e);
        mu_e = mu0*mu_rect(rect_e);
        
        for i=1:4
            grid_edge_id(i) = RectElement2Edge(rect_e,rect_i_map(i));
			grid_interface_edge_sign(i) = 1;
            if grid_edge_id(i) == 0
                grid_edge_id(i) = RectElement2AdjTriEdge(rect_e);
				grid_interface_edge_sign(i) = rect_i_sign(i);
            end
        end
        x_grid = RectElement2Grid2D(rect_e,1);
        y_grid = RectElement2Grid2D(rect_e,2);

        ExMap(x_grid,y_grid) = grid_edge_id(4);
        ExMap(x_grid,y_grid+1) = grid_edge_id(2);

        EyMap(x_grid,y_grid) = grid_edge_id(1);
        EyMap(x_grid+1,y_grid) = grid_edge_id(3);
		
		ExSign(x_grid,y_grid) = grid_interface_edge_sign(4);
		ExSign(x_grid,y_grid+1) = grid_interface_edge_sign(2);
		
		EySign(x_grid,y_grid) = grid_interface_edge_sign(1);
        EySign(x_grid+1,y_grid) = grid_interface_edge_sign(3);

        BETAx(x_grid,y_grid) = beta_x(grid_edge_id(1));
        BETAx(x_grid+1,y_grid) = beta_x(grid_edge_id(3));

        BETAy(x_grid,y_grid) = beta_y(grid_edge_id(4));
        
        BETAy(x_grid,y_grid+1) = beta_y(grid_edge_id(2));

        ALPHAx(x_grid,y_grid) = alpha_x(grid_edge_id(1));
        ALPHAx(x_grid+1,y_grid) = alpha_x(grid_edge_id(3));

        ALPHAy(x_grid,y_grid) = alpha_y(grid_edge_id(4));
        ALPHAy(x_grid,y_grid+1) = alpha_y(grid_edge_id(2));
    end
    
    if plot_logistics
        figure;
        [XGrid,YGrid] = meshgrid((1:FDTDGrid_Size-1)*DX,(1:FDTDGrid_Size-1)*DY);
        contourf(XGrid,YGrid,BETAx(2:end,:)-ALPHAx(2:end,:)+BETAy(:,2:end)-ALPHAy(:,2:end)); 
        colorbar;
        xlabel('Y (m)')
        ylabel('X (m)')
        set(gca,'fontsize',24);
    end
end
%% Assembly of the global matrix system [T] and [S].
disp('start assembling global matrix system [T] and [S] ... ')
T = zeros(Num_Edges,Num_Edges);
S = zeros(Num_Edges,Num_Edges);
f = zeros(Num_Edges,3);

ae = zeros(1,3); be = zeros(1,3); ce = zeros(1,3); 
flke = zeros(3,3);

for tri_e=1:Num_Elements_Tri
    % Element informations
    Te = zeros(3,3);Te_numerical = Te;
    Se = zeros(3,3);

    ne_tri(1:3) = Mesh_Elements_Tri(1:3,tri_e);
    where_am_I = Mesh_Nodes_Tri(:,ne_tri); 
    xe_tri = where_am_I(1,:);
    ye_tri = where_am_I(2,:);
    
    eps_e = eps0*eps_tri(tri_e);
    mu_e = mu0*mu_tri(tri_e);
    
    ae(1) = xe_tri(2)*ye_tri(3)-xe_tri(3)*ye_tri(2); ae(2) = xe_tri(3)*ye_tri(1)-xe_tri(1)*ye_tri(3); ae(3) = xe_tri(1)*ye_tri(2)-xe_tri(2)*ye_tri(1);
    be(1) = ye_tri(2) - ye_tri(3); be(2) = ye_tri(3) - ye_tri(1); be(3) = ye_tri(1) - ye_tri(2);
    ce(1) = xe_tri(3) - xe_tri(2); ce(2) = xe_tri(1) - xe_tri(3); ce(3) = xe_tri(2) - xe_tri(1);
    Delta_e = (be(1)*ce(2)-be(2)*ce(1))/2;
    
    for l=1:3
        for k=1:3
            flke(l,k)=be(l)*be(k)+ce(l)*ce(k);
        end
    end
    % assembling of local matrices
    for i=1:3 % just loop index
        edge_id_i = TriElement2Edge(tri_e,tri_i_map(i));
        edge_i_sign = TriElement2SignOfEdge(tri_e,tri_i_map(i));
        vecEdge_i = edge_i_sign*[xe_tri(tri_k_map(i))-xe_tri(tri_l_map(i)),ye_tri(tri_k_map(i))-ye_tri(tri_l_map(i))];
        edgeLength_i = norm(vecEdge_i);
        CurlOfNelk_i = edge_i_sign*CurlOfNelk2D(be(tri_l_map(i)),ce(tri_l_map(i)),be(tri_k_map(i)),ce(tri_k_map(i)),edgeLength_i,Delta_e);

        for j=1:3
            edge_id_j = TriElement2Edge(tri_e,tri_i_map(j));
            edge_j_sign = TriElement2SignOfEdge(tri_e,tri_i_map(j));
            vecEdge_j = edge_j_sign*[xe_tri(tri_k_map(j))-xe_tri(tri_l_map(j)),ye_tri(tri_k_map(j))-ye_tri(tri_l_map(j))];
            edgeLength_j = norm(vecEdge_j);
            CurlOfNelk_j = edge_j_sign*CurlOfNelk2D(be(tri_l_map(j)),ce(tri_l_map(j)),be(tri_k_map(j)),ce(tri_k_map(j)),edgeLength_j,Delta_e);
            % analytical integration for [Se] evaluation
            Se(tri_i_map(i),tri_i_map(j)) = 1/mu_e*CurlOfNelk_i*CurlOfNelk_j*Delta_e;
            
            % analytical integration for [Te] evaluation
            if tri_i_map(i)==tri_i_map(j)
                ff = edge_i_sign*edge_j_sign*eps_e*edgeLength_i*edgeLength_j/Delta_e/12;
                Te(tri_i_map(i),tri_i_map(j)) = ff*(flke(tri_l_map(i),tri_l_map(j))+flke(tri_k_map(i),tri_k_map(j))-flke(tri_l_map(i),tri_k_map(i)));
            else
                ff = edge_i_sign*edge_j_sign*eps_e*edgeLength_i*edgeLength_j/Delta_e/24;
                for k=1:3
                    if ~ismember(k,[tri_i_map(i),tri_i_map(j)])
                        break;
                    end
                end
                Te(tri_i_map(i),tri_i_map(j)) = ff*(flke(tri_l_map(i),tri_k_map(i))+flke(tri_l_map(j),tri_k_map(j))-2*flke(tri_i_map(i),tri_i_map(j))-flke(k,k));
            end
        end
    end
%     norm(Te-Te_numerical)
    % Add to global matrices
    for i=1:3
        edge_id_i = TriElement2Edge(tri_e,tri_i_map(i));
        for j=1:3
            edge_id_j = TriElement2Edge(tri_e,tri_i_map(j));
            T(edge_id_i,edge_id_j) = T(edge_id_i,edge_id_j)+Te(tri_i_map(i),tri_i_map(j)); 
            S(edge_id_i,edge_id_j) = S(edge_id_i,edge_id_j)+Se(tri_i_map(i),tri_i_map(j)); 
        end
    end
end
Conformal_RectTri_ElementPair_Count = 0;
for rect_e=1:Num_Elements_Rect
    % Element informations
    Te = zeros(4,4);
    Se = zeros(4,4);
    
    ne_rect(1:4) = Mesh_Elements_Rect(1:4,rect_e);
    where_am_I = Mesh_Nodes_Rect(:,ne_rect); 
    xe_rect = where_am_I(1,:);
    ye_rect = where_am_I(2,:);
    
    dx = xe_rect(2)-xe_rect(1);
    dy = ye_rect(4)-ye_rect(1);
    assert(dx>0&&dy>0,'Error in local node indexing: should be counter-clockwise');
    dxdy = dx*dy;
    CurlOfNelk_Rect = 1./[-dx,-dy,dx,dy];
    eps_e = eps0*eps_rect(rect_e);
    mu_e = mu0*mu_rect(rect_e);
    
    if RectElementIsAdjToBoundary(rect_e)
        tri_e = RectElement2AdjTriElement(rect_e);
        ne_tri(1:3) = Mesh_Elements_Tri(1:3,tri_e);
        where_am_I = Mesh_Nodes_Tri(:,ne_tri); 
        xe_tri = where_am_I(1,:);
        ye_tri = where_am_I(2,:);
        ae(1) = xe_tri(2)*ye_tri(3)-xe_tri(3)*ye_tri(2); ae(2) = xe_tri(3)*ye_tri(1)-xe_tri(1)*ye_tri(3); ae(3) = xe_tri(1)*ye_tri(2)-xe_tri(2)*ye_tri(1);
        be(1) = ye_tri(2) - ye_tri(3); be(2) = ye_tri(3) - ye_tri(1); be(3) = ye_tri(1) - ye_tri(2);
        ce(1) = xe_tri(3) - xe_tri(2); ce(2) = xe_tri(1) - xe_tri(3); ce(3) = xe_tri(2) - xe_tri(1);
        Delta_e = (be(1)*ce(2)-be(2)*ce(1))/2;
        
        for l=1:3
            for k=1:3
                flke(l,k)=be(l)*be(k)+ce(l)*ce(k);
            end
        end  
        
        for i=1:3
            interface_edge_id = TriElement2Edge(tri_e,tri_i_map(i));
            if EdgesIsOnBoundary(interface_edge_id)
                Conformal_RectTri_ElementPair_Count = Conformal_RectTri_ElementPair_Count + 1;
                edge_sign_tri = TriElement2SignOfEdge(tri_e,tri_i_map(i));
                edgeLength = sqrt((xe_tri(tri_k_map(i))-xe_tri(tri_l_map(i)))^2+(ye_tri(tri_k_map(i))-ye_tri(tri_l_map(i)))^2);
                CurlOfNelk_Tri = -edge_sign_tri*CurlOfNelk2D(be(tri_l_map(i)),ce(tri_l_map(i)),be(tri_k_map(i)),ce(tri_k_map(i)),edgeLength,Delta_e);
                ff = eps_e*edgeLength*edgeLength/Delta_e/12;
                temp_flke = flke(tri_l_map(i),tri_l_map(i))+flke(tri_k_map(i),tri_k_map(i))-flke(tri_l_map(i),tri_k_map(i));
                break;
            end
        end
    end

    % assembling of local matrices
    for i=1:4 % just loop index
        edge_id_i = RectElement2Edge(rect_e,rect_i_map(i));
        edge_i_sign = RectElement2SignOfEdge(rect_e,rect_i_map(i));
        vecEdge_i = edge_i_sign*[xe_rect(rect_k_map(i))-xe_rect(rect_l_map(i)),ye_rect(rect_k_map(i))-ye_rect(rect_l_map(i))];
        edgeLength_i = norm(vecEdge_i);
        CurlOfNelk_i = CurlOfNelk_Rect(rect_i_map(i));

        for j=1:4
            edge_id_j = RectElement2Edge(rect_e,rect_i_map(j));
            edge_j_sign = RectElement2SignOfEdge(rect_e,rect_i_map(j));
            vecEdge_j = edge_j_sign*[xe_rect(rect_k_map(j))-xe_rect(rect_l_map(j)),ye_rect(rect_k_map(j))-ye_rect(rect_l_map(j))];
            edgeLength_j = norm(vecEdge_j);

            CurlOfNelk_j = CurlOfNelk_Rect(rect_i_map(j));
            if edge_id_i~=0 && edge_id_j~=0 % both the test and basis funciton belong to edges of rectangular element
                % analytical integration for [Se] evaluation
                Se(rect_i_map(i),rect_i_map(j)) = 1/mu_e*CurlOfNelk_i*CurlOfNelk_j*dxdy;
                % trapezoidal integration for [Te] evaluation
                if edge_id_i==edge_id_j
                    Te(rect_i_map(i),rect_i_map(j)) = eps_e*dxdy/2;
                end
            elseif edge_id_i==0 && edge_id_j==0 % both the test and basis function belong to the edge of triangular element
                Se(rect_i_map(i),rect_i_map(j)) = 1/mu_e*CurlOfNelk_Tri*CurlOfNelk_Tri*dxdy;
                Te(rect_i_map(i),rect_i_map(j)) = ff*temp_flke*dxdy/Delta_e;
            else
                if edge_id_i==0 && edge_id_j~=0
                    Se(rect_i_map(i),rect_i_map(j)) = 1/mu_e*CurlOfNelk_Tri*CurlOfNelk_j*dxdy;
                end
                if edge_id_j==0 && edge_id_i~=0
                    Se(rect_i_map(i),rect_i_map(j)) = 1/mu_e*CurlOfNelk_Tri*CurlOfNelk_i*dxdy;
                end
            end
        end
    end
    % Add to global matrices
    for i=1:4
        edge_id_i = RectElement2Edge(rect_e,rect_i_map(i));
        for j=1:4
            edge_id_j = RectElement2Edge(rect_e,rect_i_map(j));
            if edge_id_i~=0 && edge_id_j ~=0
                T(edge_id_i,edge_id_j) = T(edge_id_i,edge_id_j)+Te(rect_i_map(i),rect_i_map(j)); 
                S(edge_id_i,edge_id_j) = S(edge_id_i,edge_id_j)+Se(rect_i_map(i),rect_i_map(j)); 
            elseif edge_id_i==0 && edge_id_j ==0
                S(interface_edge_id,interface_edge_id) = S(interface_edge_id,interface_edge_id)+Se(rect_i_map(i),rect_i_map(j)); 
            else
                if edge_id_i==0 && edge_id_j~=0
                    S(interface_edge_id,edge_id_j) = S(interface_edge_id,edge_id_j)+Se(rect_i_map(i),rect_i_map(j));
                end
                if edge_id_j==0 && edge_id_i~=0
                    S(edge_id_i,interface_edge_id) = S(edge_id_i,interface_edge_id)+Se(rect_i_map(i),rect_i_map(j));
                end
            end
        end
    end
end
assert(isequal(Conformal_RectTri_ElementPair_Count,sum(RectElementIsAdjToBoundary),sum(RectElementIsAdjToBoundary),sum(TriElementIsAdjToBoundary),sum(EdgesIsOnBoundary)),'Error in element-edge-element pairing: conformal mesh identification error')
T11 = T(1:Num_Edges_Tri,1:Num_Edges_Tri);
T22 = T(Num_Edges_Tri+1:Num_Edges,Num_Edges_Tri+1:Num_Edges);
S11 = S(1:Num_Edges_Tri,1:Num_Edges_Tri);
S22 = S(Num_Edges_Tri+1:Num_Edges,Num_Edges_Tri+1:Num_Edges);
S12 = S(1:Num_Edges_Tri,Num_Edges_Tri+1:Num_Edges);
S21 = S(Num_Edges_Tri+1:Num_Edges,1:Num_Edges_Tri);
assert(isApproximatelyEqual(S12,S21',1e-12),'Error in coupling submatrice of stiffness matrix: [S] should be symmetric.')
%% Time-stepping and Solving Linear System
disp('start time-stepping ... ')

X=linspace(-a-d1-d2,a+d1+d2,N);
Y=linspace(-a-d1-d2,a+d1+d2,N);

e_tri_coeff1 = zeros(Num_Edges_Tri,1);
e_tri_coeff2 = zeros(Num_Edges_Tri,1);
e_tri_coeff3 = zeros(Num_Edges_Tri,1);
e_rect_coeff1 = zeros(Num_Edges_Rect,1);
e_rect_coeff2 = zeros(Num_Edges_Rect,1);
e_rect_coeff3 = zeros(Num_Edges_Rect,1);

for n=1:NumTimeSteps
    t = n*dt;
    disp(['Step ',num2str(n),' at time ',num2str(t/1e-9), ' ns']);
    %% Current source in the TE plane    
    ne_tri(1:3) = Mesh_Elements_Tri(1:3,TriElementofSource);
    where_am_I = Mesh_Nodes_Tri(:,ne_tri); 
    xe_tri = where_am_I(1,:);
    ye_tri = where_am_I(2,:);
    ae(1) = xe_tri(2)*ye_tri(3)-xe_tri(3)*ye_tri(2); ae(2) = xe_tri(3)*ye_tri(1)-xe_tri(1)*ye_tri(3); ae(3) = xe_tri(1)*ye_tri(2)-xe_tri(2)*ye_tri(1);
    be(1) = ye_tri(2) - ye_tri(3); be(2) = ye_tri(3) - ye_tri(1); be(3) = ye_tri(1) - ye_tri(2);
    ce(1) = xe_tri(3) - xe_tri(2); ce(2) = xe_tri(1) - xe_tri(3); ce(3) = xe_tri(2) - xe_tri(1);
    Delta_e = (be(1)*ce(2)-be(2)*ce(1))/2;

    % evaluate the RHS of the 2nd-order updating equations
    temp_dt_list = (-1:1:1)*dt; % [previous, current, next]
    fe = zeros(1,3);
    for j = 1:3 % time step 
        temp_dt = temp_dt_list(j);
        for i = 1:3 % local edge
            jsrc_edge_id = TriElement2Edge(TriElementofSource,tri_i_map(i));
            jsrc_edge_sign = TriElement2SignOfEdge(TriElementofSource,tri_i_map(i));
            le_lk = sqrt(be(tri_i_map(i))^2+ce(tri_i_map(i))^2);
            GaussSamples = zeros(1,3);
            for ngp=1:3
                % 1. tranform the quadrature points from natural coordinates to
                % cartesian coordinates (3 points Gauss quadrature, 2nd order
                % precision) 
                x = (xe_tri(tri_l_map(ngp))+xe_tri(tri_k_map(ngp)))/2;
                y = (ye_tri(tri_l_map(ngp))+ye_tri(tri_k_map(ngp)))/2;
                % 2. evaluate functions and their dot product at the points
                Ne_i = jsrc_edge_sign*dot(Nelk2D(x,y,ae(tri_l_map(i)),be(tri_l_map(i)),ce(tri_l_map(i)),ae(tri_k_map(i)),be(tri_k_map(i)),ce(tri_k_map(i)),le_lk,Delta_e),J_dir);
                % 3. multiply with Gauss weights and scale
                GaussSamples(ngp) = -Ne_i*dJimp_dt(t+temp_dt);
            end
            fe(tri_i_map(i))=fe(tri_i_map(i))+GaussQuadTri3(GaussSamples)*2*Delta_e;
        end
        for i=1:3
            jsrc_edge_id = TriElement2Edge(TriElementofSource,tri_i_map(i));
            f(jsrc_edge_id,j)=fe(tri_i_map(i));
        end
    end
    f_tri_1=f(1:Num_Edges_Tri,1);
    f_tri_2=f(1:Num_Edges_Tri,2);
    f_tri_3=f(1:Num_Edges_Tri,3);
    
    
    %% Updating
    if use_fdtd_update_with_pml
        e_coeff2 = [e_tri_coeff2;e_rect_coeff2];
        % use conventional FDTD update in a leapfrog fashion
        % Coupling through the outer boundary edges of the truncation box
        % of FETD region
        
        Ey1 = get2DGridFrom1DCoeff(e_coeff2,EyMap,EySign);
        Ex1 = get2DGridFrom1DCoeff(e_coeff2,ExMap,ExSign);

        Hzx2 = 1 ./ BETAx(2:end,:) .* (ALPHAx(2:end,:) .* Hzx1 - mu0/(eps0*DX) * (Ey1(2:end,:)-Ey1(1:end-1,:))); % differential in x direction
        Hzy2 = 1 ./ BETAy(:,2:end) .* (ALPHAy(:,2:end) .* Hzy1 + mu0/(eps0*DY) * (Ex1(:,2:end)-Ex1(:,1:end-1))); % differential in y direction
        Hz2 =  Hzy2 + Hzx2; Hzx1 = Hzx2; Hzy1 = Hzy2;

        % enforce boundary condition of TE field: tangential H = 0
        Hz2 = [Hz2(:,1),Hz2,Hz2(:,end)]; 
        Hz2 = [Hz2(1,:);Hz2;Hz2(end,:)];
        
        Ex2 = 1 ./ BETAy .* (ALPHAy .* Ex1 + 1/DY * (Hz2(2:end-1,2:end) - Hz2(2:end-1,1:end-1))); % differential in y direction
        Ey2 = 1 ./ BETAx .* (ALPHAx .* Ey1 - 1/DX* (Hz2(2:end,2:end-1) - Hz2(1:end-1,2:end-1))); % differential in x direction
        
        e_coeff3 = update1DCoeffFrom2DGrid(e_coeff2,Ey2,EyMap);
        e_coeff3 = update1DCoeffFrom2DGrid(e_coeff3,Ex2,ExMap);
        e_rect_coeff3 = e_coeff3(Num_Edges_Tri+1:Num_Edges);
    else
        % solve a FEM linear sysetm (FDTD grid can be seen as rect element)
        % Coupling through [S21]
        K_rect = 1/dt^2*T22;
        b_rect = (2/dt^2*T22-S22)*e_rect_coeff2-1/dt^2*T22*e_rect_coeff1-S21*e_tri_coeff2;
        e_rect_coeff3 = K_rect\b_rect;
    end
    
    K_tri = (1/dt^2*T11+1/4*S11);
    b_tri = (2/dt^2*T11-1/2*S11)*e_tri_coeff2-(1/dt^2*T11+1/4*S11)*e_tri_coeff1+1/4*f_tri_3+1/2*f_tri_2+1/4*f_tri_1-S12*(1/4*e_rect_coeff3+1/2*e_rect_coeff2+1/4*e_rect_coeff1);
    % enforce PEC BC
    e_tri_coeff3(EdgesOnPEC)=0;
    e_tri_coeff3 = K_tri\b_tri;
    
    e_tri_coeff1 = e_tri_coeff2;
    e_tri_coeff2 = e_tri_coeff3;
    
    e_rect_coeff1 = e_rect_coeff2;
    e_rect_coeff2 = e_rect_coeff3;
    
    
    %% Prediction
    if apply_lpf
        t_next = t + dt;
        fe = zeros(1,3);
        for j = 1:3 % time step 
            temp_dt = temp_dt_list(j);
            for i = 1:3 % local edge
                jsrc_edge_id = TriElement2Edge(TriElementofSource,tri_i_map(i));
                jsrc_edge_sign = TriElement2SignOfEdge(TriElementofSource,tri_i_map(i));
                le_lk = sqrt(be(tri_i_map(i))^2+ce(tri_i_map(i))^2);
                GaussSamples = zeros(1,3);
                for ngp=1:3
                    % 1. tranform the quadrature points from natural coordinates to
                    % cartesian coordinates (3 points Gauss quadrature, 2nd order
                    % precision) 
                    x = (xe_tri(tri_l_map(ngp))+xe_tri(tri_k_map(ngp)))/2;
                    y = (ye_tri(tri_l_map(ngp))+ye_tri(tri_k_map(ngp)))/2;
                    % 2. evaluate functions and their dot product at the points
                    Ne_i = jsrc_edge_sign*dot(Nelk2D(x,y,ae(tri_l_map(i)),be(tri_l_map(i)),ce(tri_l_map(i)),ae(tri_k_map(i)),be(tri_k_map(i)),ce(tri_k_map(i)),le_lk,Delta_e),J_dir);
                    % 3. multiply with Gauss weights and scale
                    GaussSamples(ngp) = -Ne_i*dJimp_dt(t_next+temp_dt);
                end
                fe(tri_i_map(i))=fe(tri_i_map(i))+GaussQuadTri3(GaussSamples)*2*Delta_e;
            end
            for i=1:3
                jsrc_edge_id = TriElement2Edge(TriElementofSource,tri_i_map(i));
                f(jsrc_edge_id,j)=fe(tri_i_map(i));
            end
        end
        f_tri_1=f(1:Num_Edges_Tri,1);
        f_tri_2=f(1:Num_Edges_Tri,2);
        f_tri_3=f(1:Num_Edges_Tri,3);

        b_tri = (2/dt^2*T11-1/2*S11)*e_tri_coeff2-(1/dt^2*T11+1/4*S11)*e_tri_coeff1+1/4*f_tri_3+1/2*f_tri_2+1/4*f_tri_1-S12*(1/4*e_rect_coeff3+1/2*e_rect_coeff2+1/4*e_rect_coeff1);
        e_tri_coeff4 = K_tri\b_tri;
        e_tri_coeff2=LPF1(e_tri_coeff4,e_tri_coeff2,e_tri_coeff1);
    end
    %% Visualization
    if mod(n,plot_every_nsteps)==0
        EGrid = zeros(N,N,2);
        GridInterpolated = zeros(N,N);
        for tri_e=1:Num_Elements_Tri
            ne_tri(1:3) = Mesh_Elements_Tri(1:3,tri_e);
            where_am_I = Mesh_Nodes_Tri(:,ne_tri); 
            xe_tri = where_am_I(1,:);
            ye_tri = where_am_I(2,:);
            ae(1) = xe_tri(2)*ye_tri(3)-xe_tri(3)*ye_tri(2); ae(2) = xe_tri(3)*ye_tri(1)-xe_tri(1)*ye_tri(3); ae(3) = xe_tri(1)*ye_tri(2)-xe_tri(2)*ye_tri(1);
            be(1) = ye_tri(2) - ye_tri(3); be(2) = ye_tri(3) - ye_tri(1); be(3) = ye_tri(1) - ye_tri(2);
            ce(1) = xe_tri(3) - xe_tri(2); ce(2) = xe_tri(1) - xe_tri(3); ce(3) = xe_tri(2) - xe_tri(1);
            Delta_e = (be(1)*ce(2)-be(2)*ce(1))/2;
            for q=1:N
                for p=1:N
                    if ~GridInterpolated(q,p)
                        [in,on]=inpolygon(X(p),Y(q),xe_tri,ye_tri);
                        if (in || on) %indicating if inside the element or on the edge of the element
                            for i=1:3
                                edge_sign = TriElement2SignOfEdge(tri_e,tri_i_map(i));
                                edge_id = TriElement2Edge(tri_e,tri_i_map(i));
                                le_lk = sqrt(be(tri_i_map(i))^2+ce(tri_i_map(i))^2);
                                EGrid(q,p,:)=squeeze(EGrid(q,p,:))'+squeeze(edge_sign*Nelk2D(X(p),Y(q),ae(tri_l_map(i)),be(tri_l_map(i)),ce(tri_l_map(i)),ae(tri_k_map(i)),be(tri_k_map(i)),ce(tri_k_map(i)),le_lk,Delta_e)*e_tri_coeff2(edge_id));
                            end
                            GridInterpolated(q,p) = 1;
                        end
                    end
                end
            end
        end

        for rect_e=1:Num_Elements_Rect
            ne_rect(1:4) = Mesh_Elements_Rect(1:4,rect_e);
            where_am_I = Mesh_Nodes_Rect(:,ne_rect); 
            xe_rect = where_am_I(1,:);
            ye_rect = where_am_I(2,:);
            dx = xe_rect(2)-xe_rect(1);
            dy = ye_rect(4)-ye_rect(1);
            
            interface_edge_id = RectElement2AdjTriEdge(rect_e);
            
            C = zeros(1,4); % coeffs for local edges of rect element
            for i=1:4
                edge_id = RectElement2Edge(rect_e,rect_i_map(i));
                if edge_id~=0
                    C(rect_i_map(i)) = e_rect_coeff2(edge_id-Num_Edges_Tri);
                else
                    C(rect_i_map(i)) = e_tri_coeff2(interface_edge_id)*rect_i_sign(i);
                end
            end

            for q=1:N
                for p=1:N
                    if ~GridInterpolated(q,p)
                        [in,on]=inpolygon(X(p),Y(q),xe_rect,ye_rect);
                        if (in || on) %indicating if inside the element or on the edge of the element
                            EGrid(q,p,:)=squeeze(C(1)*N1(X(p),Y(q),xe_rect(1),dx)+C(2)*N2(X(p),Y(q),ye_rect(4),dy)+C(3)*N3(X(p),Y(q),xe_rect(2),dx)+C(4)*N4(X(p),Y(q),ye_rect(1),dy));
                            GridInterpolated(q,p) = 1;
                        end
                    end
                end
            end
        end
        figure(5);hold on;
        imagesc(X,Y,abs(real(EGrid(:,:,1)))); % x-direction on TE plane
        colormap summer;colorbar;axis image;xlabel('Y (m)');ylabel('X (m)');set(gca,'fontsize',18);
        title(['$E_x^{tot}$', ' at ', 'time = ',num2str(round(n*dt,9)),' s',],'interpreter','latex');
        for edge_id = 1:Num_Edges
            if EdgesIsOnBoundary(edge_id) 
                plot(Mesh_Nodes_Tri(1,Edge2TriNode(edge_id,:)),Mesh_Nodes_Tri(2,Edge2TriNode(edge_id,:)),'-b','linewidth',2);
            end
        end
        drawnow; 
        if save_fig
            set(gcf,'Position',[0,0,1024,1024]); 
            saveas(gca, ['Hybrid-FETD-FDTD step',num2str(n),'Ex with LPF','.png']);
        end
        pause(0.2);hold off;
    end
end
toc;

%% Functions
function [Nodes, Rectangles, NodesIsOnInnerBoundary,ElementIsInPMLRegion, ElementIsAdjToInnerBoundary, RectGrid2D] = RectangularMeshForFDTDRegion(xx,yy,R,L)
    Nx=length(xx);
    Ny=length(yy);
    NodeFlag = zeros(Nx,Ny);
    ElementInnerBoundaryFlag = zeros(Nx-1,Ny-1);
    ElementInFDTDRegionFlag = ones(Nx-1,Ny-1);
    for j=1:Ny
        for i=1:Nx
            if (abs(xx(i))+1e-12)<R && (abs(yy(j))+1e-12)<R
                NodeFlag(i,j)=nan;
            end
            if (abs(abs(xx(i))-R)<1e-12 && abs(yy(j))<=R+1e-12) || (abs(xx(i))<=R+1e-12 && abs(abs(yy(j))-R)<1e-12 )
                NodeFlag(i,j)=0.5;
            end
        end
    end
    
    e = 0;
    n = 0;
    for j=1:1:Ny-1 
        for i=1:1:Nx-1
            for jj=0:1
                for ii=0:1
                    if NodeFlag(i+ii,j+jj)==0.5
                        ElementInnerBoundaryFlag(i,j)=ElementInnerBoundaryFlag(i,j)+1;
                    end
                    if isnan(NodeFlag(i+ii,j+jj))
                        ElementInFDTDRegionFlag(i,j)=0;
                    end
                    
                end
            end
        end
    end
    % spy(ElementInFDTDRegionFlag)
    % spy(ElementInnerBoundaryFlag)
    % loop through elements
    for j=1:1:Ny-1 
        for i=1:1:Nx-1
            elementInFDTDRegionFlag = 1;
            elementInPMLRegionFlag = 0;
            
            % make nodes
            for jj=0:1
                for ii=0:1
                    if ElementInFDTDRegionFlag(i,j)
                        if NodeFlag(i+ii,j+jj)==0 || NodeFlag(i+ii,j+jj)==0.5
                            n = n+1;
                            if NodeFlag(i+ii,j+jj)==0.5
                                NodesIsOnInnerBoundary(n) = 1;
                            else
                                NodesIsOnInnerBoundary(n) = 0;
                            end
                            NodeFlag(i+ii,j+jj) = n;
                            Nodes(1,n)=xx(i+ii);
                            Nodes(2,n)=yy(j+jj);
                        end

                        if xx(i+ii)<-L && yy(j+jj)>L
                            elementInPMLRegionFlag = 1;
                        elseif abs(xx(i+ii))<L && yy(j+jj)>L
                            elementInPMLRegionFlag = 2;
                        elseif xx(i+ii)>L && yy(j+jj)>L
                            elementInPMLRegionFlag = 3;
                        elseif xx(i+ii)>L && abs(yy(j+jj))<L
                            elementInPMLRegionFlag = 4;
                        elseif xx(i+ii)>L && yy(j+jj)<-L
                            elementInPMLRegionFlag = 5;
                        elseif abs(xx(i+ii))<L && yy(j+jj)<-L
                            elementInPMLRegionFlag = 6;
                        elseif xx(i+ii)<-L && yy(j+jj)<-L
                            elementInPMLRegionFlag = 7;
                        elseif xx(i+ii)<-L && abs(yy(j+jj))<L
                            elementInPMLRegionFlag = 8;
                        end
                    end
                end
            end
            
            % make element
            if ElementInFDTDRegionFlag(i,j)
                e=e+1;
                Rectangles(1,e)=NodeFlag(i,j);
                Rectangles(2,e)=NodeFlag(i+1,j);
                Rectangles(3,e)=NodeFlag(i+1,j+1);
                Rectangles(4,e)=NodeFlag(i,j+1);
                ElementIsInPMLRegion(e)=elementInPMLRegionFlag;
                if ElementInnerBoundaryFlag(i,j)==2
                    ElementIsAdjToInnerBoundary(e)=1;
                else
                    ElementIsAdjToInnerBoundary(e)=0;
                end
                RectGrid2D(e,1) = i;
                RectGrid2D(e,2) = j;
            end
        end
    end
end
function boolean=isApproximatelyEqual(m1,m2,threshold)
    boolean = max(max(max(abs(m1-m2))))<threshold;
end
function F2=LPF1(f3,f2,F1)
    p=0.05;b=p;q=1-2*p;
    F2=b*f3+q*f2+p*F1;
end
function Grid=get2DGridFrom1DCoeff(coeff,Map,Sign)
    Grid = zeros(size(Map));
    for i=1:size(Map,1)
        for j=1:size(Map,2)
            if Map(i,j)~=0
                Grid(i,j)=coeff(Map(i,j))*Sign(i,j);
            end
        end
    end
end
function coeff2=update1DCoeffFrom2DGrid(coeff1,Grid,Map)
    coeff2 = coeff1;
    for i=1:size(Map,1)
        for j=1:size(Map,2)
            if Map(i,j)~=0
                coeff2(Map(i,j))=Grid(i,j);
            end
        end
    end
end