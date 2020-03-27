function LI = multiIndex2LinearIndex(this, MI, d, m)
%MULTIINDEX2LINEARINDEX Convert a list of multi-indices into a list of
%linear indices with lexicographical ordering
%
%   INPUT PARAMETERS:
%
%       - MI:       #Ixm list of multi-indeices
%       - d:        The degree of the mult-index (d == sum(MI))
%       - m:        The dimension of the multi-index (m == numel(MI))
%
%   OUTPUT PARAMETERS:
%
%       - LI:       #Cx1 vector of linear indices
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

% Convert combinations to linear indices ----------------------------------

LI = zeros(numInds, 1);
for i = 1:numInds
    
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

