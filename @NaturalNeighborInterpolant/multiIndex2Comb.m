function C = multiIndex2Comb(this, MI, d, m)
%MULTIINDEX2COMB Convert a list of multi-indices to a list of combinations
%
%   INPUT PARAMETERS:
%
%       - MI:       #Ixm list of multi-indeices
%       - d:        The degree of the mult-index (d == sum(MI))
%       - m:        The dimension of the multi-index (m == numel(MI))
%
%   OUTPUT PARAMETERS:
%
%       - C:        #Ixd list of combinations
%
% by Dillon Cislo 12/22/2019

% Validate Inputs ---------------------------------------------------------
validateattributes(d, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
validateattributes(m, {'numeric'}, {'scalar', 'integer', 'positive'});

validateattributes(MI, {'numeric'}, ...
    {'2d', 'ncols', m, 'nonnegative', 'integer'});

assert( all(sum(MI,2) == d), 'Multi-index has invalid degree' );

% Convert multi-indices to combinations -----------------------------------
numInds = size(MI,1); % The number of indices to convert
n = d + m - 1; % A convenience variable

C = zeros(numInds, d);
for i = 1:numInds
    
    CC = zeros(1, n);
    count = 1;
    for j = 1:m
        if MI(i,j) > 0
            CC(count:(count+MI(i,j)-1)) = 1;
            count = count+MI(i,j)+1;
        else
            count = count + 1;
        end
    end
    
   C(i,:) = find(CC);
   
end


end