function ghostPointValueHandling( this, numPoints, ghostMethod  )
%GHOSTPOINTVALUEHANDLING Calculate assigned values and gradients for the
%ghost points used in natural neighbor extrapolation.
%
%   INTEDED FOR INTERNAL USE WITH 'NaturalNeighborInterpolant.m'
%   WARNING: NO INPUT CHECKS ARE PERFORMED
%
%   INPUT PARAMETERS:
%
%       - numPoints:    The number of data points (excludes ghost points)
%
%       - ghostMethod:  The method used for determining the position of the
%                       ghost points used for natural neighbor
%                       extrapolation. The choise of method will also
%                       dictate the type of gradient estimation used
%
%   by Dillon Cislo 03/30/2020

% The vertex IDs of the ghost point
GPIDx = (numPoints+1):size(this.Points, 1);

switch ghostMethod
    
    case {'custom', 'circle', 'edge'}
        
        % Vertex IDs defining edges in the extended Delaunay triangulation
        eIDx = this.Edges;
        
        % Remove all edges containing two ghost points from the list
        GP2GP = all( ismember(eIDx, GPIDx), 2 );
        eIDx(GP2GP, :) = [];
        
        for i = 1:numel(GPIDx)
            
            % The current ghost point
            curGP = GPIDx(i);
            
            % Find any edges containing the current ghost point
            GPinE = any( eIDx == curGP, 2 );
            
            % Find the vertex IDs of the natural neighbors
            nnIDx = eIDx(GPinE, :); nnIDx(:);
            nnIDx( nnIDx == curGP ) = [];
            
            % Calculate the separation vectors between the current ghost
            % point and its natural neighbors
            X = this.Points(nnIDx,:) - ...
                repmat( this.Points(curGP,:), numel(nnIDx), 1 );
            Y = X(:,2); X = X(:,1);
            
            % Calculate the normalized inverese edge lengths
            edgeWeights = 1 ./ sqrt( X.^2 + Y.^2 );
            edgeWeights = edgeWeights ./ sum(edgeWeights);
            
            % Calculate average values and update values list
            VGP = repmat( edgeWeights, 1, size(this.Values, 2) );
            VGP = sum( VGP .* this.Values(nnIDx,:), 1 );
            this.Values(curGP,:) = VGP;
            
            % Calculate average gradients and update gradients list
            DGP = this.DataHess(nnIDx,:,1) .* X + ...
                this.DataHess(nnIDx,:,2) .* Y;
            DGP = cat(3, DGP, this.DataHess(nnIDx,:,2) .* X + ...
                this.DataHess(nnIDx,:,3) .* Y);
            DGP = DGP + this.DataGrad(nnIDx, :, :);
            DGP = repmat(edgeWeights, 1, size(this.Values, 2), 2) .* DGP;
            DGP = sum( DGP, 1 );
            this.DataGrad(curGP,:,:) = DGP;
            
            % Calculate average Hessians and update Hessian list
            HGP = repmat( edgeWeights, 1, size(this.Values, 2), 3 );
            HGP = sum( HGP .* this.DataHess(nnIDx,:,:), 1 );
            this.DataHess(curGP,:,:) = HGP;
            
        end
        
    otherwise
        
        error('Invalid ghost point construction method supplied!');
        
end

end

