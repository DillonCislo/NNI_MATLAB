function [Fq, DFq] = sibsonInterpolant(this, Xq, Yq)
%SIBSONINTERPOLANT An implementation of Sibson's C^1 continuous interpolant
%for scattered data using natural neighbor coordinates.  Output includes
%interpolated function values and analytic derivatives.  Supports
%multivariate interpolation
% 
%   INPUT PARAMETERS:
%
%       - Xq:       #Qx1 vector of the x-coordinates of query points
%
%       - Yq:       #Qx1 vector of the y-coordinates of query points
%
%   OUTPUT PARAMETERS:
%
%       - Fq:       #Qx#V matrix of interpolated function values
%
%       - DFq:      #Qx#Vx2 array of function derivatives
%
%   by Dillon Cislo 12/22/2019


%==========================================================================
% INPUT PROCESSING
%==========================================================================
% NOTE: This function is fully intended for use as an external feature.
% However, the following input checks may be commented out for speed if the
% application calls for it.

% Validate query point coordinates
validateattributes(Xq, {'numeric'}, ...
    {'vector', 'finite', 'nonnan', 'real'} );
validateattributes(Yq, {'numeric'}, ...
    {'vector', 'numel', numel(Xq), 'finite', 'nonnan', 'real'} );

if (size(Xq,2) ~= 1), Xq = Xq.'; Yq = Yq.'; end

% Check that the query points lie within the convex hull of the extended
% triangulation
[in, on] = inpolygon(Xq, Yq, this.ConvexHull(:,1), this.ConvexHull(:,2));
if any(~in)
    numOut = sum(~in);
    fmt = ['Query points: ', repmat('%g, ', 1, numOut-1), ...
        '%g lie outside the interpolation region\n'];
    fprintf(fmt, find(~in));
    error('Invalid query points supplied');
elseif any(on)
    numOn = sum(on);
    fmt = ['Query points: ', repmat('%g, ', 1, numOn-1), ...
        '%g lie on the hull of the interpolation region\n'];
    fprintf(fmt, find(~in));
    error('Invalid query points supplied');
end

% The vertices of the extended triangulation
Points = this.Points;

% The number of query points
numQueries = numel(Xq);

%==========================================================================
% CALCULATE NATURAL NEIGHBOR COORDINATES FOR EACH QUERY POINT
%==========================================================================
[allU, allUIDx, allUVC, allUA] = this.naturalNeighborCoordinates(Xq, Yq);

%==========================================================================
% CALCULATE INTERPOLATED FUNCTION VALUES AND DERIVATIVES
%==========================================================================

Fq = zeros( numQueries, size(this.Values, 2) );
DFq = zeros( numQueries, size(this.Values, 2), 2 );
for q = 1:numQueries
    
    if (this.dispInterp), progressbar(q, numQueries); end
    
    % The current query point coordinates
    vq = [Xq(q) Yq(q)];
    
    %----------------------------------------------------------------------
    % Calculate the gradients of the natural neighbor coordinates
    %----------------------------------------------------------------------
    
    % The natural neighbor coordinates of the current query point
    u = allU{q}; % #NNx1 vector
    
    % The vertex IDs of the natural neighbors of the current query point
    uIDx = allUIDx{q}; % #NNx1 vector
    
    % The coordinates of the vertices of the virtual Voronoi cell of the
    % current query point
    uVC = allUVC{q}; % #NNx2 matrix
    
    % The normalization constant of the natural neighbor coordinates
    uA = allUA(q);
    
    % Shift cell coordinates.  Used for calculated cell edge lengths
    uVC_Shift = circshift( uVC, 1, 1 ); % #NNx2 matrix
    
    % The control point function values at the natural neighbors
    ZNN = this.Values(uIDx,:); % #NNx#V matrix
    
    % The control point gradients at the natural neighbors
    GNN = this.DataGrad(uIDx,:,:); % #NNx#Vx2 array
    
    % Handle the case where a query point coincides with a control point --
    numNeighbors = numel(uIDx);
    if numNeighbors == 1
        Fq(q,:) = ZNN;
        DFq(q,:,:) = GNN;
        continue;
    end
    
    % The edge lengths of the virtual Voronoi cell. #NNx1 vector
    cellLengths = uVC - uVC_Shift;
    cellLengths = sqrt( sum( cellLengths.^2, 2 ) );
    
    % The directed edges from natural neighbor control points to the
    % query points.  NOTE: The gradient of the squared edge length is
    % simply the directed edge vector - hence the name. #NNx2 matrix
    gradL2 = repmat(vq, numNeighbors, 1) - Points(uIDx, :);
    
    % The squared/plain lengths of those edges (resp.). #NNx1 vector
    L2 = sum( gradL2.^2, 2 ); L = sqrt(L2);
    
    % The gradient of the plain edge lengths. #NNx2 matrix
    gradL = 0.5 .* gradL2 ./ L;
    
    % The gradients of the NON-NORMALIZED natural neighbor coordinates
    % with respect to the query point coordinates. A #NNx2 matrix.
    gradLU = (( uVC + uVC_Shift ) ./ 2) - repmat(vq, numNeighbors, 1);
    gradLU = cellLengths .* gradLU ./ L;
    
    % The gradients of the NORMALIZED natural neighbor coordinates
    % with respect to the query point coordinates. A #NNx2 matrix
    gradU = ( gradLU - u .* sum( gradLU, 1 ) ) ./ uA;
    
    %----------------------------------------------------------------------
    % Calculate Sibson's C^0 interpolant and its derivatives
    %----------------------------------------------------------------------
    
    % 1x#V row vector
    Z0 = sum( u .* ZNN, 1 );
    
    % 2x#V matrix
    gradZ0 = cat(1, sum(gradU(:,1) .* ZNN, 1), sum(gradU(:,2) .* ZNN, 1));
    
    %----------------------------------------------------------------------
    % Calculate the 'Beta' parameter and its derivatives
    %----------------------------------------------------------------------
    
    % A scalar
    Beta = sum( u .* L2 );
    
    % 2x1 row vector
    gradBeta = sum( gradU .* L2 + u .* gradL2, 1 ).';
    
    %----------------------------------------------------------------------
    % Calculate the 'Alpha' parameter and its derivatives
    %----------------------------------------------------------------------
    
    w = u ./ L; % #NNx1 vector
    WW = sum(w); % A scalar
    
    gradW = gradU ./ L - u .* gradL ./ L2; % #NNx2 matrix
    gradWW = sum(gradW, 1).'; % 2x1 vector
    
    % A scalar
    Alpha = sum( u .* L ) ./ WW;
    
    % 2x1 vector
    gradAlpha = sum( L .* gradU + u .* gradL, 1 ).' - Alpha .* gradWW;
    gradAlpha = gradAlpha ./ WW;
    
    %----------------------------------------------------------------------
    % Calculate the 'Xi' parameter and its derivatives
    %----------------------------------------------------------------------
    
    % #NNx#V matrix
    xi = ZNN + ( GNN(:,:,1) .* gradL2(:,1) + GNN(:,:,2) .* gradL2(:,2) );
    % NOTE: the gradient of the lower-case 'xi' = GNN --> #NNx#Vx2
    
    % 1x#V row vector
    Xi = sum( w .* xi, 1) ./ WW;

    % 2x#V matrix
    gradXi = ( permute(sum( w .* GNN, 1 ), [3 2 1]) + ...
        cat(1, sum(gradW(:,1) .* xi, 1), sum(gradW(:,2) .* xi, 1)) - ...
        gradWW .* Xi ) ./ WW;
    
    %----------------------------------------------------------------------
    % Calculate Sibson's C^1 interpolant and its derivatives
    %----------------------------------------------------------------------
    
    % A scalar
    AB = Alpha + Beta;
    
    % 1x#V row vector
    Z1 = ( Alpha .* Z0 + Beta .* Xi ) ./ AB;
    
    % 2x#V matrix
    gradZ1 = gradAlpha .* (Z0-Z1) + gradBeta .* (Xi-Z1) + ...
        Alpha .* gradZ0 + Beta .* gradXi;
    gradZ1 = gradZ1 ./ AB;
    
    %----------------------------------------------------------------------
    % Update output variables
    %----------------------------------------------------------------------
    
    Fq(q,:) = Z1;
    DFq(q,:,:) = permute(gradZ1, [3 2 1]);
        
end

end

