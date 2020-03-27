function C = linearIndex2Comb(this, LI, d, m)
%LINEARINDEX2COMB Convert a vector of linear indices with lexicographical
%ordering into a list of combinations
%
%   INPUT PARAMETERS:
%
%       - LI:       #LIx1 vector of linear indices
%       - d:        The degree of the mult-index (d == sum(MI))
%       - m:        The dimension of the multi-index (m == numel(MI))
%
%   OUTPUT PARAMETERS:
%
%       - C:        #LIxd list of combinations
%
%   by Dillon Cislo 12/22/2019

% Validate Inputs ---------------------------------------------------------
validateattributes(d, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
validateattributes(m, {'numeric'}, {'scalar', 'integer', 'positive'});

n = d + m - 1; % A convenience variable
N = nchoosek(n, d); % The maximum linear index

validateattributes(LI, {'numeric'}, ...
    {'column', 'integer', 'positive', '<=', N});

% Convert linear indices to combinations ----------------------------------
numInds = numel(LI); % The number of indices to convert

C = zeros(numInds, d);
for i = 1:numInds
    
    CC = zeros(1, d);
    
    r = 0; k = 0;
    for j = 1:(d-1)
        
        if (j ~= 1), CC(j) = CC(j-1); end
        
        while true
            
            CC(j) = CC(j)+1;
            r = nchoosek(n-CC(j), d-j);
            k = k+r;
            
            if (k >= i), break; end
            
        end
        
        k = k-r;
        
    end
    
    CC(d) = CC(d-1) + i - k;
    
    C(i,:) = CC;
    
end


end

