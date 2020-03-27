function LI = comb2LinearIndex(this, C, d, m)
%COMB2LINEARINDEX Convert a list of combinations to a list of linear
%indices sorted by lexicographical ordering
%
%   INPUT PARAMETERS:
%
%       - C:        #Cxd list of combinations
%       - d:        The degree of the mult-index (d == sum(MI))
%       - m:        The dimension of the multi-index (m == numel(MI))
%
%   OUTPUT PARAMETERS:
%
%       - LI:       #Cx1 vector of linear indices
%
%   by Dillon Cislo 12/22/2019

% Validate Inputs ---------------------------------------------------------
validateattributes(d, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
validateattributes(m, {'numeric'}, {'scalar', 'integer', 'positive'});

n = d+m-1; % A convenience variable
validateattributes(C, {'numeric'}, ...
    {'2d', 'ncols', d, '>=', 1, '<=', n, 'integer'});

% Convert combinations to linear indices ----------------------------------
numCombs = size(C,1); % The number of combinations to convert

LI = zeros(numCombs, 1);
for i = 1:numCombs
    
    Qj = 0;
    for j = 1:d
        
        if (j == 1), k1 = 1; else, k1 = C(i, j-1)+1; end
        k2 = C(i,j)-1;
        
        if k1 <= k2
            for k = k1:k2
                Qj = Qj + nchoosek(n-k, d-j);
            end
        end
        
    end

    LI(i) = 1 + Qj;
    
end



end

