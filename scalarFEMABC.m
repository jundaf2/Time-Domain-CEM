% 2d scalar (nodal) FEM formulation for scattering problem using
% first-order ABC. The finite elements are triangles with linear basis
% fucntions. Objects are PECs. Ez^tot is the primary unknown
%% settings
clear all; close all;
N = 100; % interpolation size
R = 5/3; % domain size 2R x 2R
lambda = 1; % wavelength
k0 = 2*pi/lambda; % wave number

Hmax=0.2; % maximum element size
objNum = 0; % 0, 1 or 4
obj_radius = 2/3; % < 2.5 

% the incident angle measured with respect to the positive x-axis
phi_in = 0*pi/4; 

% choose what to plot and how to plot
plot_logistics = true; % if to plot FEM matrices and meshes
use_pdeplot = false;

% Material properties and source (incident field)
eps0 = 1; % relative permittivity
mu0 = 1; % relative permeability
sigma0 = 0; % conductivity
Ez_inc_amp = 1;
Ez_inc_dir = [cos(phi_in),sin(phi_in)];

%% macros and lambda expressions
Nel2D = @(x,y,ael,bel,cel,Delta_e) 1/(2*Delta_e)*(ael+bel*x+cel*y); 
NablaNel2D = @(bel,cel,Delta_e) 1/(2*Delta_e).*[bel,cel];
GaussPoints = [-0.9061798459,-0.5384693101,0,0.5384693101,0.9061798459];
GaussWeights = [0.2369268851,0.4786286705,0.5688888889,0.4786286705,0.2369268851];
GaussQuad5 = @(GaussSamples) GaussSamples*GaussWeights';
Rotate90Mat = [cosd(90) -sind(90); sind(90) cosd(90)];
Ez_inc = @(x,y) Ez_inc_amp*exp(-1j*k0*([x,y]*Ez_inc_dir'));
%% Discretization of the 2-D domain: generate geometry & create mesh using meshpde toolbox 
model = createpde;
% physical size is measured in wavelength
%
RR = [3,4,-R,R,R,-R,-R,-R,R,R]';
% rect domain x 4 
R1 = [3,4,-R,0.5,0.5,-R,-R,-R,0.5,0.5]'; % Rectangle(square) subdomain SW
R2 = [3,4,-0.5,R,R,-0.5,-R,-R,0.5,0.5]'; % Rectangle(square) subdomain SE 
R3 = [3,4,-0.5,R,R,-0.5,-0.5,-0.5,R,R]'; % Rectangle(square) subdomain NE
R4 = [3,4,-R,0.5,0.5,-R,-0.5,-0.5,R,R]'; % Rectangle(square) subdomain NW

if objNum==4 % obj x 4 
    C1 = [1,-2.5,-2.5,1]'; % Circle object in R1
    C1 = [C1;zeros(length(R1) - length(C1),1)]; % Append extra zeros to the circles so they have the same number of rows as the rectangle
    E1 = [4,2.5,-2.5,1.5,1,0]'; % Ellipse object in R2
    E1 = [E1;zeros(length(R1) - length(E1),1)];
    r1 = [3,4,1.5,3.5,3.5,1.5,1.5,1.5,3.5,3.5]'; % Rectangle object (square) in R3
    T1 = [2,3,-3.5,-1.5,-2.5,2.5-sqrt(3)/2,2.5-sqrt(3)/2,2.5+sqrt(3)/2]'; % Triangle object in R4
    T1 = [T1;zeros(length(R1) - length(T1),1)];

    gm = [R1,R2,R3,R4,C1,E1,T1,r1];
    sf = 'R1-C1+R2-E1+R3-r1+R4-T1';
    ns = char('R1','R2','R3','R4','C1','E1','T1','r1');
elseif objNum==1
    C = [1,0,0,R/3]'; % Circle object in R1
    C = [C;zeros(length(RR) - length(C),1)]; % Append extra zeros to the circles so they have the same number of rows as the rectangle

    gm = [RR,C];
    sf = 'RR-C';
    ns = char('RR','C');
elseif objNum==0
    gm = [R1,R2,R3,R4];
    sf = 'R1+R2+R3+R4';
    ns = char('R1','R2','R3','R4');
end
ns = ns';
g = decsg(gm,sf,ns);
geometryFromEdges(model,g);
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',Hmax);
if plot_logistics == 1
    figure;
    hold on;
    pdegplot(model,'FaceLabels','on','EdgeLabels','on')
    xlim([-1.1*R,1.1*R])
    ylim([-1.1*R,1.1*R])
    xlabel('x');
    xlabel('y');
    hold off;
    
    figure;hold on;
    % pdeplot(model); 
    xlabel('x');
    xlabel('y');
    %title('Finite element discretization using first-order triangular elements')
    pdemesh(mesh,'ElementLabels','off')
    hold on
    pdemesh(mesh.Nodes,mesh.Elements,'EdgeColor','blue')
    set(gca,'fontsize',24);
    xlim([-1.1*R,1.1*R])
    ylim([-1.1*R,1.1*R])
    hold off;
end

Num_Nodes = size(mesh.Nodes,2);
Num_Elements = size(mesh.Elements,2);
X=linspace(-R,R,N);
Y=linspace(-R,R,N);

eps = ones(1,Num_Elements); % relative permittivity of the elements
mu = ones(1,Num_Elements); % relative permeability of the elements
sigma = zeros(1,Num_Elements); % conductivity

% visualization of the incident field
EzInc=zeros(N,N); 
inPEC = zeros(objNum,1);
onPEC = zeros(objNum,1);
if objNum==4 % obj x 4 
    C1 = [1,-2.5,-2.5,1]'; % Circle object in R1
    C1 = [C1;zeros(length(R1) - length(C1),1)]; % Append extra zeros to the circles so they have the same number of rows as the rectangle
    C2 = [1,2.5,-2.5,1.5]'; % Circle object in R2
    C2 = [C2;zeros(length(R1) - length(C2),1)];
    r1 = [3,4,1.5,3.5,3.5,1.5,1.5,1.5,3.5,3.5]'; % Rectangle object (square) in R3
    T1 = [2,3,-3.5,-1.5,-2.5,2.5-sqrt(3)/2,2.5-sqrt(3)/2,2.5+sqrt(3)/2]'; % Triangle object in R4
    T1 = [T1;zeros(length(R1) - length(T1),1)];
    
    for i=1:N
        for j=1:N
            [inPEC(1),onPEC(1)]=inpolygon(X(i),Y(j),[1.5,3.5,3.5,1.5],[1.5,1.5,3.5,3.5]);
            [inPEC(2),onPEC(2)]=inpolygon(X(i),Y(j),[-3.5,-1.5,-2.5],[2.5-sqrt(3)/2,2.5-sqrt(3)/2,2.5+sqrt(3)/2]);
            [inPEC(3),onPEC(3)]=incircle([X(i),Y(j)],[-2.5,-2.5],1);
            [inPEC(4),onPEC(4)]=incircle([X(i),Y(j)],[2.5,-2.5],1.5);
            if ~sum(inPEC) && ~sum(onPEC)
                EzInc(i,j) = Ez_inc(X(j),Y(i));
            end
        end
    end
elseif objNum==1
 
    for i=1:N
        for j=1:N
            [inPEC(1),onPEC(1)]=incircle([X(i),Y(j)],[0,0],R/3);

            if ~sum(inPEC) && ~sum(onPEC)
                EzInc(i,j) = Ez_inc(X(j),Y(i));
            end
        end
    end
elseif objNum==0
    for i=1:N
        for j=1:N
            if ~sum(inPEC) && ~sum(onPEC)
                EzInc(i,j) = Ez_inc(X(j),Y(i));
            end
        end
    end
end

if plot_logistics == 1
    figure;hold on;
    imagesc(X,Y,abs(real(EzInc)));
    title('E_z^{inc}')
    colormap summer
    axis image
    colorbar;
    xlabel('x');
    ylabel('y');
    set(gca,'fontsize',24);
    hold off;
    drawnow;
    pause;
end
%% Identify the source and bounaries
if objNum==4
    Nodes_onBoundary = findNodes(mesh,'region','Edge',[23,25,24,7,8,9,14,15,16,20,21,22]);
    % Find the nodes associated with edges 1, 2 and 10 (the boundary of triangle object).
    Nodes_onObj_Tri = findNodes(mesh,'region','Edge',[1,2,17]);
    % Find the nodes associated with edges 3, 4, 11 and 12 (the boundary of square object).
    Nodes_onObj_Square = findNodes(mesh,'region','Edge',[3, 4, 18, 19]);
    % Find the nodes associated with edges 20, 21, 22 and 23 (the boundary of circle object).
    Nodes_onObj_Circle = findNodes(mesh,'region','Edge',[32, 33, 34, 35]);
    % Find the nodes associated with edges 24, 25, 26, 27 and 28 (the boundary of ellipse object).
    Nodes_onObj_Ellipse = findNodes(mesh,'region','Edge',[36, 37, 38, 39, 40]);

    % for the scatter problem, PEC objectd and ABC for truncation of
    % computational domain
    Nodes_PEC = [Nodes_onObj_Tri,Nodes_onObj_Square,Nodes_onObj_Circle,Nodes_onObj_Ellipse];
    Nodes_ABC = Nodes_onBoundary;
elseif objNum==1
    Nodes_onBoundary = findNodes(mesh,'region','Edge',[1, 2, 3, 4]);

    % Find the nodes associated with edges 20, 21, 22 and 23 (the boundary of circle object).
    Nodes_onObj_Circle = findNodes(mesh,'region','Edge',[5, 6, 7, 8]);

    Nodes_PEC = Nodes_onObj_Circle;
    Nodes_ABC = Nodes_onBoundary;
elseif objNum==0
    Nodes_onBoundary = findNodes(mesh,'region','Edge',[1,2,3,7,8,9,19,20,21,22,23,24]);

    Nodes_PEC = [];
    Nodes_ABC = Nodes_onBoundary;
end


%% Assembly of the global matrix system.

K = zeros(Num_Nodes,Num_Nodes);
b = zeros(Num_Nodes,1);

kappa = 0; % curvature of the element edge --> the curvature of a straight line is zero.
gamma = kappa/2+1j*k0; % constant in mixed boundary condition, now is the parameter of first-order ABC
g = 0; % the RHS of the original PDE

for e=1:Num_Elements
    % initialization of the local element matrices and vectors
    Ke = zeros(3,3);
    Me = zeros(3,3);
    Te = zeros(3,3);
    ge = zeros(1,3);
    pe = zeros(1,3);
    xe=zeros(1,3);
    ye=zeros(1,3);
    ae=zeros(1,3);
    be=zeros(1,3);
    ce=zeros(1,3);
    
    who_am_I = mesh.Elements(1,e); % the global node idx of node 1 in element e
    where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 1 in element e
    xe(1) = where_am_I(1);
    ye(1) = where_am_I(2);
    ne(1) = who_am_I;
    who_am_I = mesh.Elements(2,e); % the global node idx of node 2 in element e
    where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 2 in element e
    xe(2) = where_am_I(1);
    ye(2) = where_am_I(2);
    ne(2) = who_am_I;
    who_am_I = mesh.Elements(3,e); % the global node idx of node 3 in element e
    where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 3 in element e
    xe(3) = where_am_I(1);
    ye(3) = where_am_I(2);
    ne(3) = who_am_I;
    
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
            Me(i,j)=-1/mu(e)*NablaNei*NablaNej'*Delta_e;
            Te(i,j)=k0^2*eps(e)*Delta_e/12*(1+(i==j));
            Ke(i,j)=Me(i,j)+Te(i,j);
        end
    end

    % Evaluation of element vector ge
    for i=1:3
        ge(i)=g*Delta_e/3;
    end
    
    % assumptions:
    % 1. the global node id number is counter clockwise
    % 2. at most one edge in an triangular element is on the boundary of
    % the computational domain
    
    % Evaluation of the element vector pe & update of the element K-matrix\
    NextLocalNodeof = [2,3,1]; % map to the next
    PreviousLocalNodeof = [3,1,2]; % map to the previous
    for lnode = 1:3 % local node
        nlnode = NextLocalNodeof(lnode); % next local node
        plnode = PreviousLocalNodeof(lnode); % previous local node

        GaussSamples = zeros(1,5);
        % if either edge 1-2, 2-3, 3-1 is on the outer boundary so that the edge is not shared by any neighboring two elements (not a common edge so that their contribution will not be canceled)
        if ismember(ne(lnode),Nodes_ABC) && ismember(ne(nlnode),Nodes_ABC) 
            vecEdge = [xe(nlnode)-xe(lnode),ye(nlnode)-ye(lnode)];
            edgeLength = norm(vecEdge);
            edgeDir = vecEdge/edgeLength;
            bn = edgeDir*Rotate90Mat; % boundary normal, rotates the array that represents the edge clockwise by 90 degrees

            % calculate the integrand value at the 5 gauss quad points 
            for ngp=1:5
                gp=GaussPoints(ngp);
                x_onEdge=(1+gp)/2*vecEdge(1)+xe(lnode);% start point, direction&&length
                y_onEdge=(1+gp)/2*vecEdge(2)+ye(lnode);
                Ne_onEdge=Nel2D(x_onEdge,y_onEdge,ae(lnode),be(lnode),ce(lnode),Delta_e);
                q_onEdge=(gamma-1j*k0*(Ez_inc_dir*bn'))*Ez_inc(x_onEdge,y_onEdge);
                GaussSamples(ngp)=-Ne_onEdge*q_onEdge;
            end
            pe(lnode)=GaussQuad5(GaussSamples)*edgeLength/2;

            % the evaluation result of the line integrals related to the
            % unknown coeffcients of the primary solution of the PDE
            for i=[lnode, nlnode]
                for j=[lnode, nlnode]
                    Ke(i,j)=Ke(i,j)-gamma*edgeLength/6*(1+(i==j));
                end
            end
        end

        if ismember(ne(lnode),Nodes_ABC) && ismember(ne(plnode),Nodes_ABC) 
            vecEdge = [xe(plnode)-xe(lnode),ye(plnode)-ye(lnode)];
            edgeLength = norm(vecEdge);
            edgeDir = vecEdge'/edgeLength;
            bn = Rotate90Mat*edgeDir; % boundary normal, rotates the array that represents the edge counter-clockwise by 90 degrees

            % calculate the integrand value at the 5 gauss quad points 
            for ngp=1:5
                gp=GaussPoints(ngp);
                x_onEdge=(1+gp)/2*vecEdge(1)+xe(lnode);% start point, direction&&length
                y_onEdge=(1+gp)/2*vecEdge(2)+ye(lnode);
                Ne_onEdge=Nel2D(x_onEdge,y_onEdge,ae(lnode),be(lnode),ce(lnode),Delta_e);
                q_onEdge=(gamma-1j*k0*(Ez_inc_dir*bn))*Ez_inc(x_onEdge,y_onEdge);
                GaussSamples(ngp)=-Ne_onEdge*q_onEdge;
            end
            pe(lnode)=GaussQuad5(GaussSamples)*edgeLength/2;

            % the evaluation result of the line integrals related to the
            % unknown coeffcients of the primary solution of the PDE
            for i=[plnode, lnode]
                for j=[plnode, lnode]
                    Ke(i,j)=Ke(i,j)-gamma*edgeLength/6*(1+(i==j));
                end
            end
        end
    end 
    
    % Add [Ke] to global [K] and [b]
    for i=1:3
        for j=1:3
            K(ne(i),ne(j)) = K(ne(i),ne(j))+Ke(i,j); 
        end
        b(ne(i))=b(ne(i))+ge(i)+pe(i);
    end
end

%% Imposition of boundary conditions.

for i=1:length(Nodes_PEC) % Dirichlet on PEC
    b(Nodes_PEC(i))=0;
    K(Nodes_PEC(i),Nodes_PEC(i))=1;
    for j=1:Num_Nodes
        if j==Nodes_PEC(i), continue; end
        b(j)=b(j)-K(j,Nodes_PEC(i))*0;
        K(Nodes_PEC(i),j)=0;
        K(j,Nodes_PEC(i))=0;
    end
end

%% Solution of the global matrix system [K]{Ez_tot}={b}.
Ez_tot=K\b;

%% Postprocessing of the results.
if use_pdeplot == true
    figure; hold on;
    % [p,e,t] = meshToPet(mesh);
    % pdeplot(p, e, t, 'XYData', abs(Ez_tot));
    pdeplot(model, 'XYData', abs(Ez_tot));
    set(gca,'fontsize',24);
    colormap summer
    axis image
    colorbar;
    title('E_z^{tot}');
    xlabel('x');
    ylabel('y');
    hold off;
else
    EzGrid = zeros(N,N);
    GridInterpolated = zeros(N,N);
    for e=1:Num_Elements
        if mod(e,100)==0
            disp(['Interpolating ',num2str(e),'th element'])
        end

        who_am_I = mesh.Elements(1,e); % the global node idx of node 1 in Elements(e)
        where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 1 in Elements(e)
        xe(1) = where_am_I(1);
        ye(1) = where_am_I(2);
        ne(1) = who_am_I;

        who_am_I = mesh.Elements(2,e); % the global node idx of node 1 in Elements(e)
        where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 1 in Elements(e)
        xe(2) = where_am_I(1);
        ye(2) = where_am_I(2);
        ne(2) = who_am_I;

        who_am_I = mesh.Elements(3,e); % the global node idx of node 1 in Elements(e)
        where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 1 in Elements(e)
        xe(3) = where_am_I(1);
        ye(3) = where_am_I(2);
        ne(3) = who_am_I;

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
                            EzGrid(q,p)=EzGrid(q,p)+Nel2D(X(p),Y(q),ae(i),be(i),ce(i),Delta_e)*Ez_tot(ne(i));
                        end
                        GridInterpolated(q,p) = 1;
                    end
                end
            end
        end
    end
    
    figure; hold on;
    title('E_z^{tot}');
    xlabel('x');
    ylabel('y');
    imagesc(X,Y,abs((EzGrid)));
    % pdemesh(mesh.Nodes,mesh.Elements);
    set(gca,'fontsize',24);
    colormap summer
    axis image
    colorbar;
    hold off;

    figure; hold on;
    title('E_z^{sca}');
    xlabel('x');
    ylabel('y');
    imagesc(X,Y,abs((EzGrid-EzInc)));
    set(gca,'fontsize',24);
    colormap summer
    axis image
    colorbar;
    hold off;

end
%% functions
function [in,on] = incircle(coord,center,radius)
    in = norm(coord-center)<radius;
    on = norm(coord-center)==radius;
end