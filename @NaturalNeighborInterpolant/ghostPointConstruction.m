function [GP, GPCH] = ghostPointConstruction( ...
    this, Xp, Yp, ghostMethod, GPin, GPe, GPr, GPn )
%GHOSTPOINTCONSTRUCTION Calculate the locations of the ghost points used
%for natural neighbor extrapolation based on the user specified
%construction method
%
%   INTEDED FOR INTERNAL USE WITH 'NaturalNeighborInterpolant.m'
%   WARNING: NO INPUT CHECKS ARE PERFORMED
%
%   INPUT PARAMETERS:
%
%       - Xp:           #Px1 list of data point x-coordinates
%
%       - Yp:           #Px2 list of data point y-coordinates
%
%       - ghostMethod:  The method used for ghost point construction
%                           - 'custom': user specified locations
%                           - 'circle': the dense bounding circle
%                           - 'edge': one ghost point per basic
%                           triangulation bounding edge
%
%       - GPin:         #GPx2 list of custom ghost point locations
%
%       - GPe:          The edge length increase factor for edge based
%                       ghost point construction
%
%       - GPr:          The radius increase factor of the ghost point
%                       circle from the circumcircle of the data point
%                       bounding box
%
%       - GPn:          The number of ghost points to create using the
%                       'circle' method
%
%   OUTPUT PARAMETERS:
%
%       - GP:           #GPx2 list of ghost point coordinates
%
%       - GPCH:         (#GP+1)x1 sorted list of ghost point IDs specifying
%                       the convex hull of the ghost points
%
%   by Dillon Cislo 03/28/2020

switch ghostMethod
    
    case 'custom'
        
        assert( ~isempty(GPin), ['Custom ghost point construction ' ...
            'selected, but no ghost points were supplied!'] );
        
        GP = GPin;
        
    case 'circle'
        
        % The bounding box of the input data points
        BBx = [ min(Xp) max(Xp) ]; BBy = [ min(Yp) max(Yp) ];
        
        % The center of the circumcircle of the bounding box
        BBCc = [ (min(Xp)+diff(BBx)/2) (min(Yp)+diff(BBy)/2) ];
        
        % The radius of the circumcircle of the bounding box
        BBCr = sqrt( diff(BBx)^2 + diff(BBy)^2 ) / 2;
        BBCr = GPr .* BBCr; % Increase radius
        
        % Calculate the locations of the ghost points
        theta = linspace( 0, 2*pi, GPn+1 )';
        theta = theta(1:(end-1));
        GP = BBCr .* [ cos(theta) sin(theta) ] + BBCc;
        
    case 'edge'
        
        % Construct basic triangulation. Note that by default boundary
        % edges of the Delaunay triangulation will be counter-clockwise
        % oriented
        TR = delaunayTriangulation( Xp, Yp );
        
        % Extract boundary edges
        fB = TR.freeBoundary;
        
        % Extract boundary edge midpoints
        Emp = [ (Xp(fB(:,2)) + Xp(fB(:,1))), ...
            (Yp(fB(:,2)) + Yp(fB(:,1))) ] ./ 2;
        
        % Extract boundary edge vectors
        Evec = [ (Xp(fB(:,2)) - Xp(fB(:,1))), ...
            (Yp(fB(:,2)) - Yp(fB(:,1))) ];
        
        % Calculate outward point edge normals
        Erot = cross( [ Evec, zeros( size(Evec,1), 1 ) ], ...
            repmat( [0 0 1], size(Evec, 1), 1), 2 );
        Erot = GPe .* Erot(:, 1:2);
        
        % Calculate the locations of the ghost points
        GP = Emp + Erot(:, 1:2);
        
    otherwise
        
        error('Invalid ghost point construction method supplied!');
end

% Check that the convex hull of the ghost points contains all
% input data points
GPCH = GP(convhull(GP), :);
if any( ~inpolygon( Xp, Yp, GPCH(:,1), GPCH(:,2) ) )
    error(['Some data points are not contained within' ...
        ' the convex hull of the ghost points']);
end


end

