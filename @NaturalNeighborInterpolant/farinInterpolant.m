function [Fq, DFq] = farinInterpolant(this, Xq, Yq)
%FARININTERPOLANT An implementation of Farins's C^1 continuous interpolant
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
%   by Dillon Cislo 12/23/2019

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

% The dimension of the interpolation function
numValues = size(this.Values, 2);

% The degree of the cubic Bezier simplex used for interpolation
d = 3;

%==========================================================================
% CALCULATE NATURAL NEIGHBOR COORDINATES FOR EACH QUERY POINT
%==========================================================================
[allU, allUIDx, allUVC, allUA] = this.naturalNeighborCoordinates(Xq, Yq);

%==========================================================================
% CALCULATE INTERPOLATED FUNCTION VALUES AND DERIVATIVES
%==========================================================================

Fq = zeros( numQueries, numValues );
DFq = zeros( numQueries, numValues, 2 );
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
    m = numel(uIDx); % The number of natural neighbors
    if m == 1
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
    gradL2 = repmat(vq, m, 1) - Points(uIDx, :);
    
    % The squared/plain lengths of those edges (resp.). #NNx1 vector
    L2 = sum( gradL2.^2, 2 ); L = sqrt(L2);
    
    % The gradients of the NON-NORMALIZED natural neighbor coordinates
    % with respect to the query point coordinates. A #NNx2 matrix.
    gradLU = (( uVC + uVC_Shift ) ./ 2) - repmat(vq, m, 1);
    gradLU = cellLengths .* gradLU ./ L;
    
    % The gradients of the NORMALIZED natural neighbor coordinates
    % with respect to the query point coordinates. A #NNx2 matrix
    gradU = ( gradLU - u .* sum( gradLU, 1 ) ) ./ uA;
    
    %----------------------------------------------------------------------
    % Calculate the coefficients of the cubic Bezier simplex
    %----------------------------------------------------------------------
    
    n = d + m - 1; % A convenience variable
    N = nchoosek(n, d); % The total number of coefficients
    
    % A matrix to store the coefficients with linear lexicographical
    % ordering
    C = zeros(N, numValues);
    
    % Calculate the linear indices of the entries of the Cij matrix -------
    LI = cell(m);
    
    for i = 1:m
        
        basisVec = zeros(1,m);
        basisVec(i) = d;
        
        LI{i,i} = basisVec;
        
        for j = 1:m
            if (i ~= j)
                
                curVec = basisVec;
                curVec(i) = curVec(i)-1;
                curVec(j) = curVec(j)+1;
                
                LI{i,j} = curVec;
                
            end
        end
       
    end
    
    % disp(cell2mat(LI(:)));
    LI = this.multiIndex2LinearIndex( cell2mat(LI(:)), d, m );
    
    % Calculate the coefficients of the Cij array -------------------------
    CIJ = zeros([m m numValues]);
    
    for i = 1:m
        
        % The function values at the current natural neighbor
        zi = ZNN(i,:); % 1x#V row vector
        
        % The function gradient values at the current natural neighbor
        gi = permute( GNN(i,:,:), [3 2 1] ); % 2x#V matrix
        
        % Control points at vertices are fixed by the interpolation
        % constraint
        CIJ(i,i,:) = zi;
        
        for j = 1:m
            if (i ~= j)
                
                % The directed edge vector from the current natural
                % neighbor to another natural neighbor. 2x1 vector
                dij = (Points(uIDx(j), :) - Points(uIDx(i), :)).';
                
                % Inner control points are constrained to be coplanar
                CIJ(i,j,:) = zi + sum(dij .* gi, 1) ./ 3;
                
            end
        end
        
    end
    
    % Update the corresponding entries of the C matrix
    C(LI, :) = reshape(CIJ, [(m*m) numValues]);
    
    % Calculate the linear indices of the Cijk coefficients ---------------
    LI = zeros( nchoosek(m,d), m );
    
    count = 0;
    for i = 1:m
        for j = (i+1):m
            for k = (j+1):m
                
                count = count+1;
                LI(count, [i j k]) = [1 1 1];
                
            end
        end
    end
    
    LI = this.multiIndex2LinearIndex( LI, d, m );
    
    % Calculate the Cijk coefficients -------------------------------------
    CIJK = zeros( nchoosek(m,d), numValues );
    
    count = 0;
    for i = 1:m
        for j = (i+1):m
            for k = (j+1):m
                
                count = count+1;
                
                uijk = ( CIJ(i,i,:) + CIJ(j,j,:) + CIJ(k,k,:) )./ 3;
                
                vijk = ( CIJ(i,j,:) + CIJ(i,k,:) + CIJ(j,i,:) + ...
                    CIJ(j,k,:) + CIJ(k,i,:) + CIJ(k,j,:) ) ./ 6;
                
                CIJK(count, :) = (3 .* vijk - uijk) ./ 2;
                
            end
        end
    end
    
    % Update the corresponding entries of the C matrix
    C(LI,:) = CIJK;
    
    %----------------------------------------------------------------------
    % Evaluate the interpolated function value by de Casteljau's method
    %----------------------------------------------------------------------
    Fq(q,:) = this.evalBezierDeCasteljau( C, u, d, m );
    
    %----------------------------------------------------------------------
    % Calculate the gradients of the interpolated function values
    %----------------------------------------------------------------------
    gradFq = zeros(m, numValues);
    
    % Generate all possible multi-indices for the gradients
    MI = nchoosek( (1:(d+m-2)), (d-1) );
    MI = this.comb2MultiIndex( MI, d-1, m );
    
    for i = 1:m
        
        % Calculate the coefficients for the current gradient -------------
        LI = MI; LI(:,i) = LI(:,i) + 1;
        LI = this.multiIndex2LinearIndex( LI, d, m );
        
        GC = C(LI,:);
        
        % Evaluate the simplex --------------------------------------------
        gradFq(i,:) = this.evalBezierDeCasteljau( GC, u, d-1, m );
        
    end
    
    gradFq = d .* gradFq;
    
    DFq(q,:,1) = sum( gradU(:,1) .* gradFq, 1 );
    DFq(q,:,2) = sum( gradU(:,2) .* gradFq, 1 );
    
end


end


% function B = evalBezierDeCasteljau( this, C, u, d, m )
% %EVALBEZIERDECASTELJAU Recursive evaluation of a Bezier-Bernstein simplex
% %by de Casteljau's algorithm.
% %
% %   INPUT PARAMETERS:
% %
% %       - C:        #Cx#V array of simplex coefficients
% %
% %       - u:        mx1 vector of barycentric coordinates
% %
% %       - d:        The degree of the simplex
% %
% %       - m:        The dimension of the multi-variate argument of the
% %                   simplex
% %
% %   OUTPUT PARAMETERS:
% %
% %       - B:        1x#V row vector of simplex function values
% 
% if size(C,1) == 1
%     
%     % BASE CASE: return the input coefficient list
%     B = C;
%     
% else
%     
%     % The number of coefficients at the next level
%     Nnew = nchoosek( d+m-2, d-1 );
%     
%     % The coefficient matrix for the next level
%     Cnew = zeros( Nnew, size(C,2) );
%     
%     % Generate all possible multi-indices for the next level
%     MI = nchoosek( (1:(d+m-2)), (d-1) );
%     MI = comb2MultiIndex( MI, d-1, m );
%     
%     for i = 1:Nnew
%         
%         % Find the linear indices of the neighbors of the new multi-indices
%         % at the previous level
%         LI = repmat(MI(i,:), m, 1) + eye(m);
%         LI = this.multiIndex2LinearIndex( LI, d, m );
%         
%         % Calculate the new coefficient 
%         Cnew(i,:) = sum( u .* C(LI,:), 1 );
%         
%     end
%     
%     B = evalBezierDeCasteljau( this, Cnew, u, d-1, m );
%     
% end
% 
% end