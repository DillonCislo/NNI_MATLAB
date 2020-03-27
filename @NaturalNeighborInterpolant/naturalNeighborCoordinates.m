function [u, uIDx, uVC, uA] = naturalNeighborCoordinates(this, Xq, Yq,...
    DelTri)
%NATURALNEIGHBORCOORDINATES Calculate Sibson's natural neighbor coordinates
%for a set of query points relative to a given set of input data points.
%It is assumed that the query poins lie within the convex hull of the data
%points. Based loosely on 'Stable Computation of Natural Neighbor
%Interpolation' by Hisamoto Hiyoshi (2008).
%
%   INTEDED FOR INTERNAL USE WITH 'NaturalNeighborInterpolant.m'
%   WARNING: NO INPUT CHECKS ARE PERFORMED
%
%   INPUT PARAMETERS:
%
%       - Xq:       #Qx1 list of the x-coordinates of query points
%
%       - Yq:       #Qx1 list of the y-coordinates of query points
%
%       - DelTri:   A Delaunay triangulation of the data points against
%                   which to calculate the natural neighbor coordinates
%
%   OUTPUT PARAMETERS:
%
%       - u:        #Qx1 cell array. Each entry contains the non-zero
%                   natural neighbor coordinates of the corresponding query
%                   point
%
%       - uIDx:     #Qx1 cell array. Each entry contains the vertex IDs of
%                   the natural neighbors of the corresponding query point
%
%       - uVC:      #Qx1 cell array. Each entry contains a
%                   counter-clockwise ordered polygon contour defining the
%                   prospective Voronoi cell of the corresponding query
%                   point
%
%       - uA:       #Qx1 vector. Each entry is the normalization constant
%                   of the natural neighbor coordinates
%
%   by Dillon Cislo 12/06/2019

%==========================================================================
% INPUT PROCESSING
%==========================================================================

if (nargin < 4)
    
    DelTri = this.DelTri;
    
    % NOTE: By default the faces of a 'delaunayTriangulation' are CCW
    % oriented
    face = this.Faces; % The vertex IDs defining faces
    vertex = this.Points; % The vertex coordinates
    
else
    
    face = DelTri.ConnectivityList;
    vertex = DelTri.Points;
    
end

numFaces = size(face, 1); % The number of faces;
numPoints = size(vertex, 1); % The number of control points

% The x-coordinates of the vertices of each face
fX = reshape( vertex(face,1), size(face) );

% The y-coordinates of the vertices of each face
fY = reshape( vertex(face,2), size(face) );

% The 'Delta' parameter of each face
if (nargin < 4)
    fDelta = this.fDelta;
else
    fDelta = this.calcDelta( fX, fY );
end

%==========================================================================
% CALCULATE NATURAL NEIGHBOR COORDINATES
%=========================================================================

u = cell( numel(Xq), 1 );
uIDx = cell( numel(Xq), 1);
uVC = cell( numel(Xq), 1);
uA = zeros( numel(Xq), 1);

for i = 1:numel(Xq)
    
    %----------------------------------------------------------------------
    % Handle the case that a query point is equal to a control point
    %----------------------------------------------------------------------
    
    % The shortest distance between the query point and any control point
    minDist = vertex - repmat([Xq(i) Yq(i)], numPoints, 1);
    minDist = sqrt(sum(minDist.^2, 2));
    
    [minDist, minIDx] = min(minDist);
    
    if (minDist < 10e-13)
        u{i} = 1; uIDx{i} = minIDx;
        continue;
    end
    
    %----------------------------------------------------------------------
    % Extract the natural neighbors of the current query point
    %----------------------------------------------------------------------
    
    % The 'Gamma' parameter of each face with respect to the query point
    fGamma = calcGamma( [fX, repmat(Xq(i), numFaces, 1)], ...
        [fY, repmat(Yq(i), numFaces, 1)] );
    
    % The 'Xi' parameter of each face with respect to the query point
    % The faces with Xi > 0 contain the query point in their circumcircle
    fXi = fGamma ./ fDelta;
    inFace = (fXi > 0);
    
    % The union of these faces is a star-shaped 'natural neighbor polygon'
    % The vertices of this polygon are the natural neighbors of the query
    % point
    nnPoly_F = face( inFace, : );
    
    % The natural neighbors of the current query point
    quIDx = unique(nnPoly_F(:));
    
    % Counter-clockwise sort the natural neighbors ------------------------
    % NOTE: It is NOT guaranteed that natural neighbor polygon is convex.
    % Sorting should be performed using original triangulation edges.
    
    % Virtual vertex IDs for the natural neighbors
    new_quIDx = (1:numel(quIDx)).';
    
    % The virtual vertex IDx defining the faces of the natural neighbor
    % polygon
    new_nnPoly_F = changem( nnPoly_F, new_quIDx, quIDx );
    
    % The free boundary of the natural neighbor polygon
    nnPoly_E = freeBoundary( triangulation( new_nnPoly_F, ...
        vertex(quIDx, :) ) );
    nnPoly_E = changem( nnPoly_E, quIDx, new_quIDx );
    
    % The sorted list of natural neighbors
    quIDx = nnPoly_E(:,1);
    
    %----------------------------------------------------------------------
    % Extract the virtual Voronoi cell of the current query point
    %----------------------------------------------------------------------
    
    % The x-coordinates of the triangles of the virtual Delaunay
    % triangulation containing the query point
    nnTri_X = reshape( vertex(nnPoly_E, 1), size(nnPoly_E) );
    nnTri_X = [ repmat(Xq(i), size(nnTri_X, 1), 1) nnTri_X ];
    
    % The x-coordinates of the triangles of the virtual Delaunay
    % triangulation containing the query point
    nnTri_Y = reshape( vertex(nnPoly_E, 2), size(nnPoly_E) );
    nnTri_Y = [ repmat(Yq(i), size(nnTri_Y, 1), 1) nnTri_Y ];
    
    % The circumcenters of those triangles are the vertices of the virtual
    % Voronoi cell of the query point
    quVC = this.triCircumcenter( nnTri_X, nnTri_Y );
    
    %----------------------------------------------------------------------
    % Calculate Sibson's natural neighbor coordinates for each natural
    % neighbor determined in the previous section
    %----------------------------------------------------------------------
    
    % A shifted list of the Voronoi vertices of the virtual cell
    quVC_Shift = circshift(quVC, 1, 1);
    
    % Calculate the circumcenters of the faces comprising the natural
    % neighbor polygons ---------------------------------------------------
    
    P_fX = fX(inFace, :); % The x-coordinates of vertices
    P_fY = fY(inFace, :); % The y-coordinates of vertices
    
    % The circumcenters of those faces
    P_CC = this.triCircumcenter( P_fX, P_fY );
    
    % Construct the coordinates -------------------------------------------    
    qA = 0; % The total area of the prospective Voronoi cell
    qu = zeros(size(quIDx)); % Vector of Sibson's coordinates
    for j = 1:numel(qu)
        
        % The circumcenters of only those faces comprising the natural
        % neighbor polygon that have the current natural neighbor as a
        % vertex
        Pj_V = P_CC( any(nnPoly_F == quIDx(j), 2), : );
        
        % Construct a counter-clockwise sorted polygon contour defining the
        % intersection between the virtual Voronoi cell of the query point
        % and the current natural neighbor. NOTE: the intersection of two
        % convex sets is ALSO convex
        Pj_V = [Pj_V; quVC(j,:); quVC_Shift(j,:)];
        Pj_V = this.sortPolygonCCW(Pj_V);
        
        % The area of the intersection region
        Pj_A = polyarea( Pj_V(:,1), Pj_V(:,2) );
        
        qA = qA + Pj_A;
        qu(j) = Pj_A;
        
    end
    
    % Normalize the coordinates
    qu = qu ./ qA;
    
    % Update output parameters
    u{i} = qu;
    uIDx{i} = quIDx;
    uVC{i} = quVC;
    uA(i) = qA;

end


end

function G = calcGamma( X, Y )
%CALCGAMMA Calculate the 'Gamma' parameter from (Hiyoshi, 2008).
%Gamma(v1, v2, v3, v4) does not have a simple geometric interpretation.
%
%   INPUT PARAMETERS:
%
%       - X:        #Nx4 list of x-coordinates
%
%       - Y:        #Nx4 list of y-coordinates
%
%   OUTPUT PARAMETERS:
%
%       - G:        #Nx1 list of 'Gamma' values

% Some convenience variables to simplify the calculation
x1 = X(:,1); x2 = X(:,2); x3 = X(:,3); x4 = X(:,4);
y1 = Y(:,1); y2 = Y(:,2); y3 = Y(:,3); y4 = Y(:,4);

% Cache some variables to avoid duplicate computations
x1S = x1.^2; x2S = x2.^2; x3S = x3.^2; x4S = x4.^2;
y1S = y1.^2; y2S = y2.^2; y3S = y3.^2; y4S = y4.^2;

x1y4y3 = x1 .* (y4 - y3);
x1y2y4 = x1 .* (y2 - y4);
x1y3y2 = x1 .* (y3 - y2);

x2y3y4 = x2 .* (y3 - y4);
x2y4y1 = x2 .* (y4 - y1);
x2y1y3 = x2 .* (y1 - y3);

x3y4y2 = x3 .* (y4 - y2);
x3y1y4 = x3 .* (y1 - y4);
x3y2y1 = x3 .* (y2 - y1);

x4y2y3 = x4 .* (y2 - y3);
x4y3y1 = x4 .* (y3 - y1);
x4y1y2 = x4 .* (y1 - y2);

% Construct the 'Gamma' vector
G = x1S .* ( x2y3y4 + x3y4y2 + x4y2y3 ) + ...
    x2S .* ( x1y4y3 + x3y1y4 + x4y3y1 ) + ...    
    x3S .* ( x1y2y4 + x2y4y1 + x4y1y2 ) + ...
    x4S .* ( x1y3y2 + x2y1y3 + x3y2y1 ) + ...
    y1S .* ( x2y3y4 + x3y4y2 + x4y2y3 ) + ...
    y2S .* ( x1y4y3 + x3y1y4 + x4y3y1 ) + ...
    y3S .* ( x1y2y4 + x2y4y1 + x4y1y2 ) + ...
    y4S .* ( x1y3y2 + x2y1y3 + x3y2y1 ) ;

end