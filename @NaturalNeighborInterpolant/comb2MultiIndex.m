function MI = comb2MultiIndex(this, C, d, m)
%COMB2MULTIINDEX Convert a list of combinations to a list of multi-indices
%
%   INPUT PARAMETERS:
%
%       - C:        #Cxd list of combinations
%       - d:        The degree of the mult-index (d == sum(MI))
%       - m:        The dimension of the multi-index (m == numel(MI))
%
%   OUTPUT PARAMETERS:
%
%       - MI:       #Cxm list of multi-indices
%
%   by Dillon Cislo 12/22/2019

% Validate Inputs ---------------------------------------------------------
validateattributes(d, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
validateattributes(m, {'numeric'}, {'scalar', 'integer', 'positive'});

n = d+m-1; % A convenience variable
validateattributes(C, {'numeric'}, ...
    {'2d', 'ncols', d, '>=', 1, '<=', n, 'integer'});

% Convert combinations to multi-indices -----------------------------------
numCombs = size(C,1); % The number of combinations to convert

MI = zeros( numCombs, m );
for i = 1:numCombs
    
    II = false(n,1);
    II(C(i,:)) = true;
    
    subs = cumsum(~II);
    subs = subs(II)+1;
    MI(i, :) = accumarray(subs, 1, [m, 1])';
    
end

end

