function [DataGrad, DataHess] = gradientEstimationDirect( this, ...
    numPoints, DVp )
%GRADIENTESTIMATIONDIRECT An implementation of a direct fitting method to
%generate gradient estimates for natural neighbor interpolation.  For a
%given set of scattered data points in the plane and a corresponding set of
%multivariate function values, gradients are chosen to minimize the
%difference between the actual function values and a 3rd order Taylor
%polynomial fit to the 3rd order natural neighborhood of each data point.
%
%   INTEDED FOR INTERNAL USE WITH 'NaturalNeighborInterpolant.m'
%   WARNING: NO INPUT CHECKS ARE PERFORMED
%
%   INPUT PARAMETERS:
%
%       - numPoints:    The number of data points (excludes ghost points)
%
%       - DVp:          #Px#Vx2 array of user specified analytic
%                       function gradients (or just an empty array)
%                       Excludes ghost point derivatives
%       
%   OUTPUT PARAMETERS:
%
%       - DataGrad:     (#P+#GP)x#Vx2 array of function gradients. Points
%                       on the convex hull will have empty entries.
%
%       - DataHess:     (#P+#GP)x#Vx3 array of function Hessians. Points
%                       on the convex hull will have empty entries.
%
% by Dillon Cislo 03/26/2020

%==========================================================================
% Input Processing
%==========================================================================

% Handle variable input arguments
if (nargin < 3), DVp = []; end

% Allocate memory for outputs
DataGrad = zeros([size(this.Values), 2]);
DataHess = zeros([size(this.Values), 3]);

% Triangulation edges
E = this.Edges;

% Remove any edges containing ghost points (their contributions will not be
% taken into account when generating gradients)
E( any(E > numPoints, 2), : ) = [];

% Generate an adjacency matrix from these edges
A = sparse( [E(:,1); E(:,2)], [E(:,2); E(:,1)], 1, numPoints, numPoints );

% Create a graph object from this adjacenct matrix
G = graph(A);

% Find distances between pairs of vertices
d = distances(G, 'Method', 'unweighted');


if isempty(DVp)
    
    %======================================================================
    % Estimate Gradients and Hessians
    %======================================================================
    
    % The inverse weights of the Taylor polynomial
    a = 1 ./ [ 1 1 2 1 2 3 2 2 3 ];
    
    for i = 1:numPoints
        
        if (this.dispGrad), progressbar(i, numPoints); end
        
        % The current point
        V = this.DelTri.Points(i,:);
        
        % Determine the 3rd order natural neighborhood of the current point
        curNN = find( (d(:,i) <= 3) & (d(:,i) > 0) );
        
        %------------------------------------------------------------------
        % Construct linear least squares problem
        %------------------------------------------------------------------
        
        % Edge vectors between neighbors and the current point
        X = this.DelTri.Points(curNN, :) - repmat(V, numel(curNN), 1);
        Y = X(:,2); X = X(:,1);
        
        % Calculate inverse edge vector lengths
        invL = 1 ./ sqrt( X.^2 + Y.^2 );
        
        % The least-squares LHS matrix
        A = repmat( a, numel(curNN), 1 ) .* repmat( invL, 1, 9 );
        A = A .* [ X, Y, X.^2, X .* Y, Y.^2, ...
            X.^3, X.^2 .* Y, X .* Y.^2, Y.^3 ];
        
        %------------------------------------------------------------------
        % Solve the linear system for the gradient of each value component
        %------------------------------------------------------------------
        
        for j = 1:size(this.Values, 2)
            
            % The least-squares RHS vector
            B = invL .* ( this.Values(curNN, j) - ...
                repmat( this.Values(i,j), numel(curNN), 1 ) );
            
            % Calculate the gradient for the current component
            g = A \ B;
            
            % Update the output gradient array
            DataGrad(i,j,1) = g(1);
            DataGrad(i,j,2) = g(2);
            
            % Update the output Hessian array
            DataHess(i,j,1) = g(3);
            DataHess(i,j,2) = g(4);
            DataHess(i,j,3) = g(5);
            
        end
        
    end
    
else
    
    %======================================================================
    % Estimate Hessians From Analytic Gradients
    %======================================================================
    
    % Set analytic derivatives
    DataGrad(1:numPoints, :, :) = DVp;
    
    % The inverse weights of the truncated Taylor polynomial
    a = 1 ./ [ 2 1 2 3 2 2 3 ];
    
    for i = 1:numPoints
        
        if (this.dispGrad), progressbar(i, numPoints); end
        
        % The current point
        V = this.DelTri.Points(i,:);
        
        % Determine the 3rd order natural neighborhood of the current point
        curNN = find( (d(:,i) <= 3) & (d(:,i) > 0) );
        
        %------------------------------------------------------------------
        % Construct linear least squares problem
        %------------------------------------------------------------------
        
        % Edge vectors between neighbors and the current point
        X = this.DelTri.Points(curNN, :) - repmat(V, numel(curNN), 1);
        Y = X(:,2); X = X(:,1);
        
        % Calculate inverse edge vector lengths
        invL = 1 ./ sqrt( X.^2 + Y.^2 );
        
        % The least-squares LHS matrix
        A = repmat( a, numel(curNN), 1 ) .* repmat( invL, 1, 7 );
        A = A .* [ X.^2, X .* Y, Y.^2, ...
            X.^3, X.^2 .* Y, X .* Y.^2, Y.^3 ];
        
        %------------------------------------------------------------------
        % Solve the linear system for the gradient of each value component
        %------------------------------------------------------------------
        
        for j = 1:size(this.Values, 2)
            
            taylorPol = this.Values(i,j) + ...
                DataGrad(i,j,1) .* X + DataGrad(i,j,2) .* Y;
            % taylorPol = repmat( taylorPol, numel(curNN), 1 );
            
            % The least-squares RHS vector
            B = invL .* ( this.Values(curNN, j) - taylorPol );
            
            % Calculate the gradient for the current component
            g = A \ B;
            
            % Update the output Hessian array
            DataHess(i,j,1) = g(1);
            DataHess(i,j,2) = g(2);
            DataHess(i,j,3) = g(3);
            
        end
        
    end
    
    
end

end

