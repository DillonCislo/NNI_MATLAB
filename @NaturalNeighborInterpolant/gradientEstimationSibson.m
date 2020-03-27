function DataGrad = gradientEstimationSibson( this, numPoints )
% GRADIENTESTIMATIONSIBSON An implementation of Sibson's C^1 continuous
% gradient estimation method.  For a given set of scattered data points in
% the plane and a corresponding set of multivariate function values,
% gradients are estimated for points within the convex hull from the least
% squares plane through values at each points natural neighbors
%
%   INTEDED FOR INTERNAL USE WITH 'NaturalNeighborInterpolant.m'
%   WARNING: NO INPUT CHECKS ARE PERFORMED
%
%   INPUT PARAMETERS:
%
%       - numPoints:    The number of data points (excludes ghost points)
%
%   OUTPUT PARAMETERS:
%
%       - DataGrad:     (#P+#GP)x#Vx2 array of function gradients. Points
%                       on the convex hull will have empty entries.
%
% by Dillon Cislo 12/17/2019

DataGrad = zeros([size(this.Values), 2]);

% An anonymous function used to map indices in the reduced triangulation to
% the original triangulation
indexMap = @(i,j) i + (i > j);

for i = 1:numPoints
    
    if (this.dispGrad), progressbar(i, numPoints); end
    
    % The original triangulation (includes ghost points)
    DT = this.DelTri;
    
    % The coordinates of the current data point
    V = DT.Points(i,:);
    
    % The reduced triangulation without the current data point
    DT.Points(i,:) = [];
    
    %----------------------------------------------------------------------
    % Calculate the natural neighbor coordinates of the current data point
    %----------------------------------------------------------------------
    [u, uIDx, ~, ~] = this.naturalNeighborCoordinates( V(1), V(2), DT );
    u = u{1}; uIDx = uIDx{1};
    
    %----------------------------------------------------------------------
    % Calculate the least-squares weights and matrix
    %----------------------------------------------------------------------
    
    % The directed edge vectors
    eIJ = DT.Points( uIDx, : ) - repmat( V, numel(u), 1 );
    
    % The squared edge vector lenghts
    lIJ2 = sum( eIJ.^2, 2 );
    
    % The least-squares weights
    wIJ = sqrt( u ./ lIJ2 );
    
    % The least-squares LHS matrix
    A = wIJ .* eIJ;
    
    %----------------------------------------------------------------------
    % Solve the linear system for the gradient of each value component
    %----------------------------------------------------------------------
    
    % The indices of the natural neighbors in the original triangulation
    uIDx0 = indexMap(uIDx, i);
    
    for j = 1:size(this.Values, 2)
        
        % The least-squares RHS vector ------------------------------------
        B = wIJ .* ( this.Values(uIDx0, j) - ...
            repmat(this.Values(i,j), numel(u), 1) );
        
        % Calculate the gradient for the current component
        g = A \ B;
        
        % Update the output array
        DataGrad(i,j,1) = g(1);
        DataGrad(i,j,2) = g(2);
        
    end
    
end
    
end
