% VECTOR FINITE ELEMENT METHOD
% 
clear all; close all; clc; tic;
%% settings
N = 100; % interpolation size
Hmax=0.1; % element size
plot_logistics = false; % if to plot FEM matrices and meshes

lambda=1; 
k0 = 2*pi/lambda;

a = lambda*2/3; % radius of circle object
d = lambda*1; % distance between obj and truncation box

% the incident angle measured with respect to the positive x-axis
phi_in = 0*pi/4; 
E_inc_amp = 1;
E_inc_dir = [cos(phi_in),sin(phi_in), 0];
E_field_dir = [0,1,0]/sqrt(1);
inspectDir = [0,1,0]/sqrt(1);
%% mappings, macros and lambda expressions
i_map = [3,1,2]; % from loop index to local edge id
l_map = [1,2,3]; % from loop index to local start node id
k_map = [2,3,1]; % from loop index to local end node id
Rotate90Mat = [cosd(90) -sind(90); sind(90) cosd(90)];

GaussWeightsTri3 = [1/3,1/3,1/3];
GaussQuadTri3 = @(GaussSamples) GaussSamples*GaussWeightsTri3';

GaussPointsLine5 = [-0.9061798459,-0.5384693101,0,0.5384693101,0.9061798459];
GaussWeightsLine5 = [0.2369268851,0.4786286705,0.5688888889,0.4786286705,0.2369268851];
GaussQuadLine5 = @(GaussSamples) GaussSamples*GaussWeightsLine5';

Nel2D = @(x,y,ael,bel,cel,Delta_e) 1/(2*Delta_e)*(ael+bel*x+cel*y); 
NablaNel2D = @(ael,bel,cel,Delta_e) 1/(2*Delta_e).*[bel,cel];
Nelk2D = @(x,y,ael,bel,cel,aek,bek,cek,le_lk,Delta_e) (Nel2D(x,y,ael,bel,cel,Delta_e)*NablaNel2D(aek,bek,cek,Delta_e)-Nel2D(x,y,aek,bek,cek,Delta_e)*NablaNel2D(ael,bel,cel,Delta_e))*le_lk;
CurlOfNelk2D = @(bel,cel,bek,cek,le_lk,Delta_e) le_lk/(2*Delta_e^2)*(bel*cek-bek*cel);

E_inc = @(x,y) E_inc_amp*exp(-1j*k0*([x,y,0]*E_inc_dir'))*E_field_dir; % vector field
CurlOfEinc = @(x,y) E_inc_amp*exp(-1j*k0*([x,y,0]*E_inc_dir'))*1j*k0*cross(E_field_dir,E_inc_dir); % \hat{E} \cross \arrow{k0} curl of 3D vector field

%% visualization of a vector finite element basis and its curl in 2D
% using a circle of radius 1 in polar coordinate
if plot_logistics
    figure;hold on;
    % Generate mesh in x \in [-1,1], y \in [-1,1]
    x_range = -1:0.05:1;
    y_range = x_range;
    [X,Y] = meshgrid(x_range,y_range);
    
    % information of the finite element
    scaling_factor = 0;
    while scaling_factor<0.1*0.6
        rand_theta = sort(2*pi*rand(1,3));
        rand_x = 1*cos(rand_theta);
        rand_y = 1*sin(rand_theta);
        xe1 = rand_x(1);xe2 = rand_x(2);xe3 = rand_x(3);
        ye1 = rand_y(1);ye2 = rand_y(2);ye3 = rand_y(3);
        ae1 = xe2*ye3-xe3*ye2;
        ae2 = xe3*ye1-xe1*ye3;
        ae3 = xe1*ye2-xe2*ye1;
        be1 = ye2 - ye3;% y 2-3
        be2 = ye3 - ye1;% y 3-1
        be3 = ye1 - ye2;% y 1-2
        ce1 = xe3 - xe2;% x 3-2
        ce2 = xe1 - xe3;% x 1-3
        ce3 = xe2 - xe1;% x 2-1
        Delta_e = (be1*ce2-be2*ce1)/2;
        Nabla_Nel = NablaNel2D(ae1,be1,ce1,Delta_e);
        Nabla_Nek = NablaNel2D(ae2,be2,ce2,Delta_e);
        le_12 = sqrt(be3^2+ce3^2);
        scaling_factor = 0.1*min([sqrt(be1^2+ce1^2),sqrt(be2^2+ce2^2),sqrt(be3^2+ce3^2)]);
    end
    % Curl of Basis
    CurlOfNelk = CurlOfNelk2D(be1,ce1,be2,ce2,le_12,Delta_e);
    AbsReCur = abs(CurlOfNelk);
    IntAbsReCur = AbsReCur / (double(vpa(1,2))/(1+floor(1 ./ 10.^floor(log10(1)))));
    max_color_value = 100*ceil(max(IntAbsReCur));
    jet_color = colormap(summer(ceil(max_color_value)));
    selected_color = jet_color(ceil(100*IntAbsReCur),:);
    patch(rand_x,rand_y,selected_color);
    
    for p=1:length(x_range)
        for q=1:length(y_range)
            [in,on]=inpolygon(X(p,q),Y(p,q),rand_x,rand_y);
            in_or_on = in||on;
            if in_or_on
                % Basis itself
                Ne_l = Nel2D(X(p,q),Y(p,q),ae1,be1,ce1,Delta_e);
                Ne_k = Nel2D(X(p,q),Y(p,q),ae2,be2,ce2,Delta_e);
                vXvY = (Ne_l*Nabla_Nek-Ne_k*Nabla_Nel)*le_12;
                
                vX=vXvY(1);vY=vXvY(2);
                QUIV = quiver(X(p,q),Y(p,q),vX*scaling_factor,vY*scaling_factor,'k','linewidth',2);
                QUIV.ShowArrowHead = 'off';
                QUIV.Marker = '.';
                
                % Divergence of Basis
            end
        end
    end
    
    axis([-1 1 -1 1],'square');
    grid on
    colormap gray;
    hold off;
end

% sum of the three
xe = zeros(1,3);
ye = zeros(1,3);
ae = zeros(1,3);
be = zeros(1,3);
ce = zeros(1,3);
if plot_logistics
    figure;hold on;
    % Generate mesh in x \in [-1,1], y \in [-1,1]
    x_range = -1:0.05:1;
    y_range = x_range;
    [X,Y] = meshgrid(x_range,y_range);
    
    % information of the finite element
    scaling_factor = 0;
    while scaling_factor<0.1*0.6
        rand_theta = sort(2*pi*rand(1,3));
        rand_x = 1*cos(rand_theta);
        rand_y = 1*sin(rand_theta);
        xe(1) = rand_x(1);xe(2) = rand_x(2);xe(3) = rand_x(3);
        ye(1) = rand_y(1);ye(2) = rand_y(2);ye(3) = rand_y(3);
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
        scaling_factor = 0.1*min([sqrt(be(1)^2+ce(1)^2),sqrt(be(2)^2+ce(2)^2),sqrt(be(3)^2+ce(3)^2)]);
    end
    patch(rand_x,rand_y,[1,1,1]);
    for p=1:length(x_range)
        for q=1:length(y_range)
            [in,on]=inpolygon(X(p,q),Y(p,q),rand_x,rand_y);
            in_or_on = in||on;
            if in_or_on
                % sum-up three basis functions
				vXvY = [0,0]; % init
                for i=1:3
                    Nabla_Nel = NablaNel2D(ae(l_map(i)),be(l_map(i)),ce(l_map(i)),Delta_e);
                    Nabla_Nek = NablaNel2D(ae(k_map(i)),be(k_map(i)),ce(k_map(i)),Delta_e);
                    le_lk = sqrt(be(i_map(i))^2+ce(i_map(i))^2);
                    Ne_l = Nel2D(X(p,q),Y(p,q),ae(l_map(i)),be(l_map(i)),ce(l_map(i)),Delta_e);
                    Ne_k = Nel2D(X(p,q),Y(p,q),ae(k_map(i)),be(k_map(i)),ce(k_map(i)),Delta_e);
                    vXvY = vXvY + (Ne_l*Nabla_Nek-Ne_k*Nabla_Nel)*le_lk;
                end
				vX=vXvY(1);vY=vXvY(2);
				QUIV = quiver(X(p,q),Y(p,q),vX*scaling_factor,vY*scaling_factor,'k','linewidth',2);
				QUIV.ShowArrowHead = 'off';
				QUIV.Marker = '.';
            end
        end
    end
    
    axis([-1 1 -1 1],'square');
    grid on
    colormap gray;
    hold off;
end

%% Visualization of the incident field
X=linspace(-a-d,a+d,N);
Y=linspace(-a-d,a+d,N);
EIncGrid=zeros(N,N); 

for i=1:N
    for j=1:N
        [inPEC,onPEC]=incircle([X(i),Y(j)],[0,0],a);
        if ~inPEC && ~onPEC
            Evec = E_inc(X(j),Y(i));
            EIncGrid(i,j) = dot(Evec,inspectDir); % project to a arbitrary dimension to get scalar field
        end
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
[p,e,t] = meshToPet(mesh);

if plot_logistics == 1
    figure; hold on;
    pdegplot(g,'FaceLabels','on','EdgeLabels','on')
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
Num_Elements = size(t,2); % number of elements�

eps = ones(1,Num_Elements); % relative permittivity of the elements
mu = ones(1,Num_Elements); % relative permeability of the elements
sigma = zeros(1,Num_Elements); % conductivity

%% Edge numbering and its connection to nodes and elements
% loop through the elements and give each edge its connected nodes
% edge to node, the first node is the start of an edge while the second
% is the end of an edge
% Element-to-Edge Connectivity Array for a Triangular Mesh
% 2D Array: row\column --> Num_Elements*3 e \ ne(1, 2; e) ne(1, 3; e) ne(2,
% 3; e) --> ne(l, k; e)
Num_Edges = 0;
Element2Edge = zeros(Num_Elements,3);
Edge2Node = zeros(3*Num_Elements,2); % [global_edge_id,global_node_id], do not have a specific direction
Element2SignOfEdge = zeros(Num_Elements,3);
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
%% Assembly of the global matrix system [T], [R] and [S].
K = zeros(Num_Edges,Num_Edges); % global LHS matrix
b = zeros(Num_Edges,1);   % global RHS vector
for e=1:Num_Elements
    % initialization of the local element matrices and vectors
    Me = zeros(3,3);
    Te = zeros(3,3);
    Se = zeros(3,3);
    pe = zeros(1,3);
    
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
    eta_er = sqrt(mu(e)/eps(e));
    
    elementIsOnBoundary = 0; % flag to tell if this element is on boundary
    for i=1:3 % just loop index
        edge_i_sign = Element2SignOfEdge(e,i_map(i));
        edge_id_i = Element2Edge(e,i_map(i)); % global edge id for local edge i_map(i)
        vecEdge_i = [xe(k_map(i))-xe(l_map(i)),ye(k_map(i))-ye(l_map(i))];
        edgeLength_i = norm(vecEdge_i);
        edgeDir_i = vecEdge_i/edgeLength_i;
        bni2D = edgeDir_i*Rotate90Mat; % boundary normal, rotates the array that represents the edge clockwise by 90 degrees
        bni3D = [bni2D,0];
        CurlOfNelk_i = edge_i_sign*CurlOfNelk2D(be(l_map(i)),ce(l_map(i)),be(k_map(i)),ce(k_map(i)),edgeLength_i,Delta_e);

        for j=1:3
            edge_j_sign = Element2SignOfEdge(e,i_map(j));
            edge_id_j = Element2Edge(e,i_map(j)); % global edge id for edge j 
            vecEdge_j = [xe(k_map(j))-xe(l_map(j)),ye(k_map(j))-ye(l_map(j))];
            edgeLength_j = norm(vecEdge_j);
            edgeDir_j = vecEdge_j/edgeLength_j;
            bnj2D = edgeDir_j*Rotate90Mat; % boundary normal, rotates the array that represents the edge clockwise by 90 degrees
            bnj3D = [bnj2D,0];
            CurlOfNelk_j = edge_j_sign*CurlOfNelk2D(be(l_map(j)),ce(l_map(j)),be(k_map(j)),ce(k_map(j)),edgeLength_j,Delta_e);
            % analytical integration for [Me] evaluation
            Me(i_map(i),i_map(j)) = CurlOfNelk_i*CurlOfNelk_j*Delta_e;
            
            % 2D numerical integration for [Te] evaluation
            GaussSamples = zeros(1,3);
            for ngp=1:3
                % 1. tranform the quadrature points from natural coordinates to
                % cartesian coordinates (3 points Gauss quadrature, 2nd order
                % precision) 
                x = (xe(l_map(ngp))+xe(k_map(ngp)))/2;
                y = (ye(l_map(ngp))+ye(k_map(ngp)))/2;
                % 2. evaluate functions and their dot product at the points
                Nelk_i = edge_i_sign*Nelk2D(x,y,ae(l_map(i)),be(l_map(i)),ce(l_map(i)),ae(k_map(i)),be(k_map(i)),ce(k_map(i)),edgeLength_i,Delta_e);
                Nelk_j = edge_j_sign*Nelk2D(x,y,ae(l_map(j)),be(l_map(j)),ce(l_map(j)),ae(k_map(j)),be(k_map(j)),ce(k_map(j)),edgeLength_j,Delta_e);
                % 3. multiply with Gauss weights and scale
                GaussSamples(ngp) = dot(Nelk_i,Nelk_j);
            end
            Te(i_map(i),i_map(j)) = GaussQuadTri3(GaussSamples)*2*Delta_e;
            
            if EdgesIsOnABC(edge_id_i) && EdgesIsOnABC(edge_id_j) && (i==j) % i and j are the same and on ABC
                elementIsOnBoundary = 1;
                % 1D numerical integration for [Se] evaluation
                GaussSamples = zeros(1,5);
                for ngp=1:5
                    gp=GaussPointsLine5(ngp);
                    x_onEdge=(1+gp)/2*vecEdge_j(1)+xe(l_map(j));% start point, direction&&length
                    y_onEdge=(1+gp)/2*vecEdge_j(2)+ye(l_map(j));
                    Nelk_i_onEdge2D = edge_i_sign*Nelk2D(x_onEdge,y_onEdge,ae(l_map(i)),be(l_map(i)),ce(l_map(i)),ae(k_map(i)),be(k_map(i)),ce(k_map(i)),edgeLength_i,Delta_e);
                    Nelk_i_onEdge3D = [Nelk_i_onEdge2D,0];
                    Nelk_j_onEdge2D = edge_j_sign*Nelk2D(x_onEdge,y_onEdge,ae(l_map(j)),be(l_map(j)),ce(l_map(j)),ae(k_map(j)),be(k_map(j)),ce(k_map(j)),edgeLength_j,Delta_e);
                    Nelk_j_onEdge3D = [Nelk_j_onEdge2D,0];
                    GaussSamples(ngp)=dot(cross(bnj3D,Nelk_i_onEdge3D),cross(bnj3D,Nelk_j_onEdge3D));
                end
                Se(i_map(i),i_map(j))=GaussQuadLine5(GaussSamples)*edgeLength_j/2;
            end
            if EdgesIsOnABC(edge_id_j) && (i~=j) % edge j is on ABC
                elementIsOnBoundary = 1;
                % 1D numerical integration for [Se] evaluation
                GaussSamples = zeros(1,5);
                for ngp=1:5
                    gp=GaussPointsLine5(ngp);
                    x_onEdge=(1+gp)/2*vecEdge_j(1)+xe(l_map(j));% start point, direction&&length
                    y_onEdge=(1+gp)/2*vecEdge_j(2)+ye(l_map(j));
                    Nelk_i_onEdge2D = edge_i_sign*Nelk2D(x_onEdge,y_onEdge,ae(l_map(i)),be(l_map(i)),ce(l_map(i)),ae(k_map(i)),be(k_map(i)),ce(k_map(i)),edgeLength_i,Delta_e);
                    Nelk_i_onEdge3D = [Nelk_i_onEdge2D,0];
                    Nelk_j_onEdge2D = edge_j_sign*Nelk2D(x_onEdge,y_onEdge,ae(l_map(j)),be(l_map(j)),ce(l_map(j)),ae(k_map(j)),be(k_map(j)),ce(k_map(j)),edgeLength_j,Delta_e);
                    Nelk_j_onEdge3D = [Nelk_j_onEdge2D,0];
                    GaussSamples(ngp)=dot(cross(bnj3D,Nelk_i_onEdge3D),cross(bnj3D,Nelk_j_onEdge3D));
                end
                Se(i_map(i),i_map(j))=GaussQuadLine5(GaussSamples)*edgeLength_j/2;
            end
            if EdgesIsOnABC(edge_id_i) && (i~=j) % edge i is on ABC
                elementIsOnBoundary = 1;
                % 1D numerical integration for [Se] evaluation
                GaussSamples = zeros(1,5);
                for ngp=1:5
                    gp=GaussPointsLine5(ngp);
                    x_onEdge=(1+gp)/2*vecEdge_i(1)+xe(l_map(i));% start point, direction&&length
                    y_onEdge=(1+gp)/2*vecEdge_i(2)+ye(l_map(i));
                    Nelk_i_onEdge2D = edge_i_sign*Nelk2D(x_onEdge,y_onEdge,ae(l_map(i)),be(l_map(i)),ce(l_map(i)),ae(k_map(i)),be(k_map(i)),ce(k_map(i)),edgeLength_i,Delta_e);
                    Nelk_i_onEdge3D = [Nelk_i_onEdge2D,0];
                    Nelk_j_onEdge2D = edge_j_sign*Nelk2D(x_onEdge,y_onEdge,ae(l_map(j)),be(l_map(j)),ce(l_map(j)),ae(k_map(j)),be(k_map(j)),ce(k_map(j)),edgeLength_j,Delta_e);
                    Nelk_j_onEdge3D = [Nelk_j_onEdge2D,0];
                    GaussSamples(ngp)=dot(cross(bni3D,Nelk_i_onEdge3D),cross(bni3D,Nelk_j_onEdge3D));
                end
                Se(i_map(i),i_map(j))=GaussQuadLine5(GaussSamples)*edgeLength_i/2;
            end
        end
        
        % numerical integration for vector [pe] evaluation
        if EdgesIsOnABC(edge_id_i) %(EdgesIsOnABC(edge_id_i) || (ismember(Edge2Node(edge_id_i,1),NodesOnBoundary) || ismember(Edge2Node(edge_id_i,2),NodesOnBoundary))) && elementIsOnBoundary            
            GaussSamples = zeros(1,5);
            % calculate the integrand value at the 5 gauss quad points 
            for ngp=1:5
                gp=GaussPointsLine5(ngp);
                x_onEdge=(1+gp)/2*vecEdge_i(1)+xe(l_map(i));% start point, direction&&length
                y_onEdge=(1+gp)/2*vecEdge_i(2)+ye(l_map(i));
                KN_onEdge = cross(bni3D,1/mu(e)*CurlOfEinc(x_onEdge,y_onEdge))+1j*k0/eta_er*cross(bni3D,cross(bni3D,E_inc(x_onEdge,y_onEdge)));
                Nelk_onEdge2D = edge_i_sign*Nelk2D(x_onEdge,y_onEdge,ae(l_map(i)),be(l_map(i)),ce(l_map(i)),ae(k_map(i)),be(k_map(i)),ce(k_map(i)),edgeLength_i,Delta_e);
                Nelk_onEdge3D = [Nelk_onEdge2D,0];
                GaussSamples(ngp)=-dot(Nelk_onEdge3D,KN_onEdge);
            end
            pe(i_map(i))=pe(i_map(i))+GaussQuadLine5(GaussSamples)*edgeLength_i/2;
        end
    end

    Ke=1/mu(e)*Me-k0^2*eps(e)*Te+1j*k0/eta_er*Se;
    
    % Add [Ke] and [pe] to global [K] and [b]
    for i=1:3
        edge_id_i = Element2Edge(e,i_map(i));
        for j=1:3
            edge_id_j = Element2Edge(e,i_map(j));
            K(edge_id_i,edge_id_j) = K(edge_id_i,edge_id_j)+Ke(i_map(i),i_map(j)); 
        end
        b(edge_id_i)=b(edge_id_i)+pe(i_map(i));
    end
end

%% Imposition of PEC boundary conditions.
for edge_id=1:Num_Edges % Dirichlet on PEC
    if EdgesIsOnPEC(edge_id)
        b(edge_id)=0;
        K(edge_id,:)=0;
        K(:,edge_id)=0;
        K(edge_id,edge_id)=1;
    end
end
%% Solve Linear system
Etot_coefficients = K\b;
%% Post-processing
EtotGrid = zeros(N,N); 
Etot = zeros(N,N,2); % basis are vectors
GridInterpolated = zeros(N,N);
for e=1:Num_Elements
    if mod(e,100)==0
        disp(['Interpolating ',num2str(e),'th element'])
    end
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
                x = X(p);
                y = Y(q);
                [in,on]=inpolygon(x,y,xe,ye);
                if (in || on) 
                    for i=1:3
                        edge_sign = Element2SignOfEdge(e,i_map(i));
                        edge_id = Element2Edge(e,i_map(i));
                        le_lk = sqrt(be(i_map(i))^2+ce(i_map(i))^2);
                        Etot(q,p,:)=squeeze(Etot(q,p,:))'+squeeze(edge_sign*Nelk2D(x,y,ae(l_map(i)),be(l_map(i)),ce(l_map(i)),ae(k_map(i)),be(k_map(i)),ce(k_map(i)),le_lk,Delta_e)*Etot_coefficients(edge_id));
                    end
                    GridInterpolated(q,p) = 1;
                end
            end
        end
    end
end

for i=1:N
    for j=1:N
        Evec = [squeeze(Etot(i,j,:))',0];
        EtotGrid(i,j) = dot(Evec,inspectDir); % project to a arbitrary dimension to get scalar field
    end
end

%% Visualization
figure; hold on;
title('E^{tot}');
xlabel('x');
ylabel('y');
imagesc(X,Y,abs((EtotGrid)),'Interpolation','bilinear');
set(gca,'fontsize',24);
colormap summer
axis image
colorbar;
hold off;

figure; hold on;
title('E^{sca}');
xlabel('x');
ylabel('y');
imagesc(X,Y,abs((EtotGrid-EIncGrid)),'Interpolation','bilinear');
set(gca,'fontsize',24);
colormap summer
axis image
colorbar;
hold off;

figure;hold on;
imagesc(X,Y,abs(real(EIncGrid)),'Interpolation','bilinear');
title('E^{inc}')
colormap summer
axis image
colorbar;
xlabel('x');
ylabel('y');
set(gca,'fontsize',24);
hold off;
toc;
%% functions
function [in,on] = incircle(coord,center,radius)
    in = norm(coord-center)<radius;
    on = norm(coord-center)==radius;
end