function B = evalBezierDeCasteljau( this, C, u, d, m )
%EVALBEZIERDECASTELJAU Recursive evaluation of a Bezier-Bernstein simplex
%by de Casteljau's algorithm.
%
%   INPUT PARAMETERS:
%
%       - C:        #Cx#V array of simplex coefficients
%
%       - u:        mx1 vector of barycentric coordinates
%
%       - d:        The degree of the simplex
%
%       - m:        The dimension of the multi-variate argument of the
%                   simplex
%
%   OUTPUT PARAMETERS:
%
%       - B:        1x#V row vector of simplex function values

if size(C,1) == 1
    
    % BASE CASE: return the input coefficient list
    B = C;
    
else
    
    % The number of coefficients at the next level
    Nnew = nchoosek( d+m-2, d-1 );
    
    % The coefficient matrix for the next level
    Cnew = zeros( Nnew, size(C,2) );
    
    % Generate all possible multi-indices for the next level
    MI = nchoosek( (1:(d+m-2)), (d-1) );
    MI = this.comb2MultiIndex( MI, d-1, m );
    
    for i = 1:Nnew
        
        % Find the linear indices of the neighbors of the new multi-indices
        % at the previous level
        LI = repmat(MI(i,:), m, 1) + eye(m);
        LI = this.multiIndex2LinearIndex( LI, d, m );
        
        % Calculate the new coefficient 
        Cnew(i,:) = sum( u .* C(LI,:), 1 );
        
    end
    
    B = this.evalBezierDeCasteljau( Cnew, u, d-1, m );
    
end

end
