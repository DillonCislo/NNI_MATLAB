function [DataGrad, DataHess] = gradientEstimationSibson( this, ...
    numPoints, iterAlpha, DVp )
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
%       - iterAlpha:    A scalar on [0,1] that control the relative
%                       importance of each subproblem in the second stage
%                       of iterative derivative estimation.  Values close
%                       to 0 try to match second iteration derivatives to
%                       first iteration derivatives.  Values close to 1 try
%                       to match second iteration derivatives to function
%                       values
%
%       - DVp:          #Px#Vx2 array of user specified analytic
%                       function gradients (or just an empty array).
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
% by Dillon Cislo 12/17/2019

% Handle variable input arguments
if (nargin < 4), DVp = []; end

% Allocate memory for outputs/temporary storage
DataGrad1 = zeros([size(this.Values), 2]);
DataGrad = zeros([size(this.Values), 2]);
DataHess = zeros([size(this.Values), 3]);

% An anonymous function used to map indices in the reduced triangulation to
% the original triangulation
indexMap = @(i,j) i + (i >= j);

%==========================================================================
% Calculate Natural Neighbor Coordinates in Reduced Triangulations
%==========================================================================
u = cell(numPoints, 1);
uIDx = cell(numPoints, 1);
uA = zeros(numPoints, 1);

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
    [curU, curUIDx, ~, curUA] = ...
        this.naturalNeighborCoordinates( V(1), V(2), DT );
    
    curU = curU{1};
    curUIDx = indexMap(curUIDx{1}, i);
    curUA = curUA(1);
    
    % Prune ghost points and renormalize ----------------------------------
    ghostNeighbors = curUIDx > numPoints;
    
    if any( ghostNeighbors )
        
        curUIDx(ghostNeighbors) = [];
        
        curU = curUA .* curU;
        curU(ghostNeighbors) = [];
        
        curUA = sum(curU);
        curU = curU ./ curUA;
        
    end
    
    % Update global cell arrays
    u{i} = curU;
    uIDx{i} = curUIDx;
    uA(i) = curUA;
    
end

if isempty( DVp )
    
    %======================================================================
    % Estimate Gradients and Hessians
    %======================================================================
    
    %----------------------------------------------------------------------
    % STEP 1: Estimate gradients usings Sibson's method
    %----------------------------------------------------------------------
    
    for i = 1:numPoints
        
        % Construct least-squares weights ---------------------------------
        
        % The location of the current point
        V = this.DelTri.Points(i,:);
        
        % The directed edge vectors between neighbors and the current point
        X = this.DelTri.Points( uIDx{i}, : ) - ...
            repmat( V, numel(u{i}), 1 );
        
        Y = X(:,2);
        X = X(:,1);
        
        % The squared edge vector lenghts
        lIJ2 = X.^2 + Y.^2;
        
        % The least-squares weights
        wIJ = sqrt( u{i} ./ lIJ2 );        
        
        % Construct least-squares problem ---------------------------------
        
        % The least-squares LHS matrix
        A = wIJ .* [X Y];
        
        % Solve least-squares problem for each value component ------------
        
        for j = 1:size(this.Values, 2)
            
            % The least-squares RHS vector
            B = wIJ .* ( this.Values(uIDx{i}, j) - ...
                repmat(this.Values(i,j), numel(u{i}), 1) );
            
            % Calculate the gradient for the current component
            g = A \ B;
            
            % Update the output  gradient array
            DataGrad1(i,j,1) = g(1);
            DataGrad1(i,j,2) = g(2);
            
        end
        
    end
    
    %----------------------------------------------------------------------
    % STEP 2: Estimate gradients and Hessians using iterative fit
    %----------------------------------------------------------------------
    
    for i = 1:numPoints
        
        % Construct least-squares weights ---------------------------------
        
        % The location of the current point
        V = this.DelTri.Points(i,:);
        
        % The directed edge vectors between neighbors and the current point
        X = this.DelTri.Points( uIDx{i}, : ) - ...
            repmat( V, numel(u{i}), 1 );
        
        Y = X(:,2);
        X = X(:,1);
        
        % The squared edge vector lenghts
        lIJ2 = X.^2 + Y.^2;
        
        % The least-squares weights
        wIJ = sqrt( u{i} ./ lIJ2 );
        
        % Construct the LHS matrix of each sub-problem --------------------
        
        Z = zeros(size(X));
        I = ones(size(X));
        
        % The least-squares LHS matrix of the first sub-problem
        A1 = iterAlpha .* repmat( wIJ, 1, 5 );
        A1 = A1 .* [ X, Y, (X.^2 ./ 2), (X .* Y), (Y.^2 ./ 2) ];
        
        % The least-squares LHS matrix of the second sub-problem
        A2 = (1-iterAlpha) .* repmat(wIJ, 2, 5);
        A2 = A2 .* [ I Z X Y Z; Z I Z X Y ];
        
        % Combine to form the full least-squares LHS matrix
        A = [A1; A2];
        
        % Solve the complete least-squares problem for each component -----
        
        for j = 1:size(this.Values, 2)
            
            % The least-squares RHS vector of the first sub-problem
            B1 = iterAlpha .* wIJ .* ( this.Values(uIDx{i}, j) - ...
                repmat(this.Values(i,j), numel(u{i}), 1) );
            
            % The least-squares RHS vector of the second sub-problem
            B2 = (1-iterAlpha) .* repmat(wIJ, 2, 1);
            B2 = B2 .* [ DataGrad1(uIDx{i},j,1); DataGrad1(uIDx{i},j,2) ];
            
            % Combine to form the full least-squares RHS vector
            B = [B1; B2];
            
            % Calculate the gradients/Hessian for the current component
            g = A \ B;
            
            % Update output gradient array
            DataGrad(i,j,1) = g(1);
            DataGrad(i,j,2) = g(2);
            
            % Update output Hessian array
            DataHess(i,j,1) = g(3);
            DataHess(i,j,2) = g(4);
            DataHess(i,j,3) = g(5);
            
        end
        
    end
    
else
    
    %======================================================================
    % Estimate Hessians From Analytic Gradients
    %======================================================================
    % WARNING: THIS METHOD IS NOT RECOMMENDED AS THERE IS NO GUARANTEE (OR
    % EVEN HEURISTIC SUGGESTION) THAT THE LEAST-SQUARES PROBLEM WILL NOT BE
    % RANK DEFICIENT
    %
    % If analytic gradients are supplied without analytic Hessian data it
    % is recommended that the user instead use the
    % 'gradientEstimationDirect' functionality
    
    % Set analytic derivates
    DataGrad(1:numPoints, :, :) = DVp;
    
    for i = 1:numPoints
        
        % Construct least-squares weights ---------------------------------
        
        % The location of the current point
        V = this.DelTri.Points(i, :);
        
        % The directed edge vectors between neighbors and the current point
        X = this.DelTri.Points( uIDx{i}, : ) - ...
            repmat( V, numel(u{i}), 1 );
        
        Y = X(:,2);
        X = X(:,1);
        
        % The squared edge vector lenghts
        lIJ2 = X.^2 + Y.^2;
        
        % The least-squares weights
        wIJ = sqrt( u{i} ./ lIJ2 );
        
        % Construct least-squares problem ---------------------------------
        
        % The least-squares LHS matrix
        A = repmat(wIJ, 1, 3) .* [ (X.^2 ./ 2), (X .* Y), (Y.^2 ./ 2) ];
        
        % Solve least squares problem for each value component ------------
        
        for j = 1:size(this.Values, 2)
            
            taylorPol = this.Values(i,j) + ...
                DataGrad(i,j,1) .* X + DataGrad(i,j,2) .* Y;
            
            % The least-squares RHS vector
            B = wIJ .* ( this.Values(uIDx{i}, j) - taylorPol );
            
            % Calculate the Hessian for the current component
            g = A \ B;
            
            % Update output Hessian array
            DataHess(i,j,1) = g(1);
            DataHess(i,j,2) = g(2);
            DataHess(i,j,3) = g(3);
            
        end
        
    end
    
end

end
