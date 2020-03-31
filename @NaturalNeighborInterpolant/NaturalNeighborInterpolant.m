classdef NaturalNeighborInterpolant < handle
    %NATURALNEIGHBORINTERPOLANT An implementation of Sisbon's and Farin's
    %C^1 continuous interpolants for scattered data using natural neighbor
    %coordinates.  Supports C^1 continuous extrapolation using the
    %assigned-value Ghost Point method. Output includes interpolated values
    %with analytic derivates (if you don't want derivatives just used
    %MATLAB's 'scatteredInterpolant!).  Supports vector valued
    %interpolation.
    %
    % by Dillon Cislo 12/04/2019
    
    %======================================================================
    %======================================================================
    %                           PROPERTIES
    %======================================================================
    %======================================================================
    
    properties (SetAccess = protected)
        
        % (#P+#GP)x#V matrix.  Holds the function values of the
        % scattered data and ghost points.
        Values
        
        % (#P+#GP)x2 matrix.  Holds the scattered data and ghost point
        % locations
        Points
        
        % #Fx3 matrix.  Holds the connectivity list of the extended
        % Delaunay triangulation
        Faces
        
        % #Ex2 matrix. Holds the edge connectivity list of the extended
        % Delaunay triangulation
        Edges
        
        % The extended Delaunay triangulation of the scattered 
        % data and ghost points
        DelTri
        
        % (#P+#GP)x#Vx2 array. Holds the gradient data of the
        % input and ghost points
        DataGrad
        
        % (#P+#GP)x#Vx3 array. Holds the Hessian data of the
        % input and ghost points
        DataHess
        
        % The 'Delta' parameter of each face of the extended Delaunay
        % triangulation.  Equal to twice the signed area of the face
        fDelta
        
        % #GPx2 matrix.  Holds the coordinates of the convex hull of the
        % extended Delaunay triangulaiton
        ConvexHull
        
        % Display progress towards estimating discrete gradients
        dispGrad
        
        % Display progress towards interpolating values
        dispInterp
        
        
    end
    
    %======================================================================
    %======================================================================
    %                              METHODS
    %======================================================================
    %======================================================================
    
    methods (Access = public)
        
        function this = NaturalNeighborInterpolant(Xp, Yp, Vp, varargin)
            %NATURALNEIGHBORINTERPOLANT Default constructor
            %
            %   MANDATORY INPUT PARAMETERS:
            %
            %       - Xp:       #Px1 list of data point x-coordinates
            %
            %       - Yp:       #Px2 list of data point y-coordinates
            %
            %       - Vp:       #Px#V matrix of data point function values
            %
            %   OPTIONAL INPUT PARAMETERS ((Name, Value)-Pairs):
            %
            %       - {'GhostPointMethod', ghostMethod}: The method used
            %       for determining the position of the ghost points used
            %       for natural neighbor extrapolation. The choice of
            %       method will also dictacte the type of gradient
            %       estimation used.
            %           - ghostMethod:  'circle'
            %
            %       - {'GhostPoints', GP}: Allows the user to explicitly
            %       specify the chost points used for extrapolation.
            %       Default is the dense circumcircle of the bounding box
            %       of the data point with a radius increase factor x2
            %           - GP:   #GPx2 list of ghost point coordinates
            %
            %       - {'GhostPointEdgeFactor', GPe}: The edge length
            %       increase factor for edge based ghost point construction
            %           - GPe: {1}
            %
            %       - {'GhostPointRadusFactor', GPr}: The radius increase
            %       factor of the ghost point circle from the circumcircle
            %       of the data point bounding box
            %           - GPr: {2}
            %
            %       - {'GhostPointNumber', GPn}: The number of ghost points
            %           - GPn: { min(10*numel(Xp), 5000) }
            %
            %       - {'Gradients', DVp}: The analytic gradients of the
            %       data point function
            %           - DVp: #Px#Vx2 array
            %
            %       - {'Hessians', HVp}: The analytic Hessians of the data
            %       point function
            %           - HVp: #Px#Vx3 array
            %
            %       - {'GradientEstimation', gradType}: The method to use
            %       for gradient estimation
            %           - gradType: 'Direct'
            %
            %       - {'GradIterAlpha', iterAlpha}: The relative importance
            %       of the two least-squares subproblems in the iterative
            %       form of Sibson's method for gradient estimation
            %           - iterAlpha: 0.001
            %
            %       - {'DisplayGradProgress', dispGrad}: If true
            %       discrete gradient estimation progress will be displayed
            %           - dispGrad: {true}
            %
            %       - {'DisplayInterpProgress', dispInterp}; if true
            %       interpolation progress will be displayed
            %           - dispInterp: {true}
            %
            %   OUTPUT PARAMETERS:
            %
            %       - this:      An instance of this class
            
            %--------------------------------------------------------------
            % Input Processing
            %--------------------------------------------------------------
            
            % Validate Mandatory Inputs -----------------------------------
            validateattributes(Xp, {'numeric'}, ...
                {'vector', 'finite', 'nonnan', 'real'} );
            validateattributes(Yp, {'numeric'}, ...
                {'vector', 'numel', numel(Xp), ...
                'finite', 'nonnan', 'real'} );
            validateattributes(Vp, {'numeric'}, ...
                {'2d', 'nrows', numel(Xp), 'finite', 'nonnan', 'real'} );
            
            if (size(Xp,2) ~= 1), Xp = Xp.'; Yp = Yp.'; end
            
            % Validate Optional Inputs ------------------------------------
            
            % Default option values
            GP = [];
            ghostMethod = 'circle';
            GPe = 1;
            GPr = 2;
            GPn = min( 10*numel(Xp), 2000 );
            DVp = [];
            HVp = [];
            gradientSupplied = false;
            hessianSupplied = false;
            gradType = 'direct';
            iterAlpha = 0.001;
            dispGrad = true;
            dispInterp = true;
            
            for i = 1:length(varargin)
                if isa(varargin{i},'double')
                    continue;
                end
                if isa(varargin{i},'logical')
                    continue;
                end
                if ~isempty(regexp(varargin{i}, ...
                        '^[Gg]host[Pp]oints', 'match'))
                    GP = varargin{i+1};
                    validateattributes(GP, {'numeric'}, ...
                        {'2d', 'ncols', 2, 'finite', 'nonnan', 'real'});
                end
                if ~isempty(regexp(varargin{i}, ...
                        '^[Gg]host[Pp]oint[Mm]ethod', 'match'))
                    ghostMethod = lower(varargin{i+1});
                    validateattributes(ghostMethod, {'char'}, ...
                        {'scalartext'});
                end
                if ~isempty(regexp(varargin{i}, ...
                        '^[Gg]host[Pp]oint[Ee]dge[Ff]actor', 'match'))
                    GPe = varargin{i+1};
                    validateattributes(GPe, {'numeric'}, ...
                        {'scalar', 'positive', 'real', 'finite'});
                end
                if ~isempty(regexp(varargin{i}, ...
                        '^[Gg]host[Pp]oint[Rr]adius[Ff]actor', 'match'))
                    GPr = varargin{i+1};
                    validateattributes(GPr, {'numeric'}, ...
                        {'scalar', 'positive', 'real', 'finite'});
                    if (GPr <= 1)
                        error(['Ghost point radius increase factor' ...
                            ' must be greater than one']);
                    end
                end
                if ~isempty(regexp(varargin{i}, ...
                        '^[Gg]host[Pp]oint[Nn]umber', 'match'))
                    GPn = varargin{i+1};
                    validateattributes(GPn, {'numeric'}, ...
                        {'scalar', 'positive', 'integer', ...
                        'real', 'finite'});
                end
                if ~isempty(regexp(varargin{i}, ...
                        '^[Gg]radients', 'match'))
                    DVp = varargin{i+1};
                    if ~isequal(size(DVp), [numel(Xp) size(Vp,2) 2])
                        error('Gradient input is improperly sized');
                    end
                    gradientSupplied = true;
                end
                if ~isempty(regexp(varargin{i}, ...
                        '^[Hh]essians', 'match'))
                    HVp = varargin{i+1};
                    if ~isequal(size(HVp), [numel(Xp) size(Vp,2) 3])
                        error('Hessian input is improperly sized');
                    end
                    hessianSupplied = true;
                end
                if ~isempty(regexp(varargin{i}, ...
                        '^[Gg]radient[Ee]stimation', 'match'))
                    gradType = lower(varargin{i+1});
                    validateattributes(gradType, {'char'}, {'scalartext'});
                end
                if ~isempty(regexp(varargin{i}, ...
                        '^[Gg]rad[Ii]ter[Aa]lpha', 'match'))
                    iterAlpha = varargin{i+1};
                    validateattributes(iterAlpha, {'numeric'}, ...
                        {'scalar', 'real', 'finite', '<=', 1, '>=' 0});
                end
                if ~isempty(regexp(varargin{i}, ...
                        '^[Dd]isplay[Gg]rad[Pp]rogress', 'match'))
                    dispGrad = varargin{i+1};
                    validateattributes(dispGrad, ...
                        {'logical'}, {'scalar'});
                end
                if ~isempty(regexp(varargin{i}, ...
                        '^[Dd]isplay[Ii]nterp[Pp]rogress', 'match'))
                    dispInterp = varargin{i+1};
                    validateattributes(dispInterp, ...
                        {'logical'}, {'scalar'});
                end
                    
            end
            
            % Display output handling
            this.dispGrad = dispGrad;
            this.dispInterp = dispInterp;
            
            %--------------------------------------------------------------
            % Ghost Point Construction
            %--------------------------------------------------------------
            
            if ~isempty(GP)
                ghostMethod = 'custom';
                GPn = size(GP,1);
            end
            
            [GP, ~] = this.ghostPointConstruction( Xp, Yp, ...
                ghostMethod, GP, GPe, GPr, GPn );
            
            GPn = size(GP,1);
            
            this.Points = [ Xp Yp; GP ];
            
            %--------------------------------------------------------------
            % Build Delaunay Triangulation/Calculate 'Delta' Parameter
            %--------------------------------------------------------------
            
            % Build triangulation
            this.DelTri = delaunayTriangulation( this.Points );
            this.Faces = this.DelTri.ConnectivityList;
            this.Edges = this.DelTri.edges;
            
            % Determine convex hull
            DTCH = this.DelTri.convexHull;
            DTCH(end, :) = [];
            
            % Check that the convex hull only contains ghost points
            assert( all( DTCH > numel(Xp) ), ...
                'Extended convex hull contains true data points!');
            
            this.ConvexHull = this.Points(DTCH(1:(end-1)), :);
            
            % The x-coordinates of the vertices of each face
            fX = reshape( this.Points(this.Faces,1), size(this.Faces) );
            
            % The y-coordinates of the vertices of each face
            fY = reshape( this.Points(this.Faces,2), size(this.Faces) );
            
            this.fDelta = this.calcDelta( fX, fY );
            
            %--------------------------------------------------------------
            % Gradient Estimation
            %--------------------------------------------------------------
            
            % The function values for non-ghost points must be set prior to
            % gradient estimation
            this.Values = [ Vp; zeros([GPn size(Vp,2)]) ];
            
            if gradientSupplied
                
                if hessianSupplied
                    
                    % Set analytic derivatives
                    DataGrad = cat( 1, DVp, zeros([GPn size(Vp,2) 2]) );
                    DataHess = cat( 1, HVp, zeros([GPn size(Vp,2) 2]) );
                    
                else
                    
                    % Estimate Hessians for all points within the
                    % convex hull
                    if dispGrad, fprintf('Estimating gradients: '); end
                    
                    if strcmp( gradType, 'direct' )
                        
                        [DataGrad, DataHess] = ...
                            this.gradientEstimationDirect( ...
                            numel(Xp), DVp );
                        
                    elseif strcmp( gradType, 'sibson' )
                        
                        [DataGrad, DataHess] = ...
                            this.gradientEstimationSibson( ...
                            numel(Xp), iterAlpha, DVp );
                        
                    else
                        
                        error(['Invalid gradient estimation ' ...
                            'method supplied!']);
                        
                    end
                    
                end
                
            else
                
                if hessianSupplied
                    
                    error(['Gradient estimation from analytic ' ...
                        'Hessian data is not supported!']);
                    
                else
                    
                    % Estimate Hessians for all points within the
                    % convex hull
                    if dispGrad, fprintf('Estimating gradients: '); end
                    
                    if strcmp( gradType, 'direct' )
                        
                        [DataGrad, DataHess] = ...
                            this.gradientEstimationDirect(numel(Xp));
                        
                    elseif strcmp( gradType, 'sibson' )
                        
                        [DataGrad, DataHess] = ...
                            this.gradientEstimationSibson( ...
                            numel(Xp), iterAlpha );
                        
                    else
                        
                        error(['Invalid gradient estimation ' ...
                            'method supplied!']);
                        
                    end
                    
                end
                
            end
            
            % Set gradients/Hessian for non-ghost points
            this.DataGrad = DataGrad;
            this.DataHess = DataHess;
            
            % Validate gradients/Hessian
            validateattributes( this.DataGrad, {'numeric'}, ...
                {'finite', 'nonnan', 'real'} );
            validateattributes( this.DataHess, {'numeric'}, ...
                {'finite', 'nonnan', 'real'} );
                        
            %--------------------------------------------------------------
            % Ghost Point Value Handling
            %--------------------------------------------------------------
            % Values assigned to ghost points are just the weighted average
            % of the values at the non-ghost point natural neighbors.
            % Weights are the normalized inverse distances to the natural
            % neighbors.  If analytic gradients were supplied gradient
            % values are averaged similarly.
            
            this.ghostPointValueHandling( numel(Xp), ghostMethod );
            
        end
        
        function [ ccwPoly, order ] = sortPolygonCCW(this, poly)
            %SORTPOLYGONCCW Sort the vertices of a convex planar polygon
            %in counter-clockwise order
            %
            %   INPUT PARAMETERS:
            %
            %       - poly:     #Px2 list of unordered polygon vertices
            %
            %   OUTPUT PARAMETERS:
            %
            %       - ccwPoly:  #Px2 list of sorted polygon vertices
            %
            %       - order:    #Px1 list of indices into the old ordering
            %                   that produces the new ordering
            
            % The center of mass of the polygon verices
            COM = mean(poly, 1);
            
            % Shift the polgon vertices so that the center of mass
            % is at the origin
            polyShift = poly - COM;
            
            % Find the angles of the shifted polygon vertices
            ang = atan2( polyShift(:,2), polyShift(:,1) );
            
            % Find the correct sorted order
            [~, order] = sort(ang);
            
            % Re order the polygon vertices
            ccwPoly = poly(order, :);
            
        end
        
        function CC = triCircumcenter(this, X, Y )
            %CIRCUMCENTER A vectorized calculation of the circumcenters of
            %a set of planar triangles
            %
            %   INPUT PARAMETERS:
            %
            %       - X:        #Tx3 list of triangle vertex x-coordinates
            %
            %       - Y:        #Tx3 list of triangle vertex y-coordinates
            %
            %   OUTPUT PARAMETERS:
            %
            %       - CC:       #Tx2 list of circumcenter coordinates
            
            % Some convenience variables to simplify the calculation
            x1 = X(:,1); x2 = X(:,2); x3 = X(:,3);
            y1 = Y(:,1); y2 = Y(:,2); y3 = Y(:,3);
            
            % Caching variables to avoid duplicate computations
            y12   = y2-y1; y31   = y1-y3; y23   = y3-y2;
            y12_P = y2+y1; y31_P = y1+y3; y23_P = y3+y2;
            y123 = y12 .* y23 .* y31;
            
            x1S = x1.^2; x2S = x2.^2; x3S = x3.^2;
            x12S = x2S - x1S;
            x23S = x3S - x2S;
            x31S = x1S - x3S;
            
            x3y12 = x3 .* y12;
            x2y31 = x2 .* y31;
            x1y23 = x1 .* y23;
            
            % Calculate circumcenter coordinates --------------------------
            
            % The numerator for the x-coordinates
            xNum = x3 .* x3y12 + x2 .* x2y31 + x1 .* x1y23 - y123;
            
            % The numerator for the y-coordinates
            yNum = x3 .* x12S + x2 .* x31S + x1 .* x23S + ...
                x3y12 .* y12_P + x2y31 .* y31_P + x1y23 .* y23_P;
                
            % The denominator for both coordinates
            denom = 2.* ( x3y12 + x2y31 + x1y23 );
            
            CC = [ xNum yNum ] ./ denom;
            
        end

        
    end
    
    methods (Access = protected)
        
        function D = calcDelta( this, X, Y )
            %CALCDELTA Calculate the 'Delta' parameter from 
            %(Hiyoshi, 2008). Delta(v1, v2, v3) is equal to twice the 
            %signed area of the triangle defined by the points {v1, v2, v3}
            %
            %   INPUT PARAMETERS:
            %
            %       - X:        #Nx3 list of x-coordinates
            %
            %       - Y:        #Nx3 list of y-coordinates
            %
            %   OUTPUT PARAMETERS:
            %
            %       - D:        #Nx1 list of 'Delta' values
            
            D = X(:,3) .* ( Y(:,1) - Y(:,2) ) + ...
                X(:,1) .* ( Y(:,2) - Y(:,3) ) + ...
                X(:,2) .* ( Y(:,3) - Y(:,1) );
            
        end
        
    end
    
    methods (Access = private)
        
        
    end
    
end

