%**************************************************************************
% Draw random fractions (fraction) out of a given distribution (MAT)
%**************************************************************************
function MAT = subset(MAT,fraction)
%==========================================================================

% Call:   MAT = subset(MAT,fraction)
%         example:  subset(MAT,0.5) --> draws a random 50% subset out of
%         distribution MAT
%
% Input:  MAT:          distribution, column vector
%         fraction:     0<fraction<1, fraction of randomly drawn subsample
%
% Output: MAT:          subsample
%	
%   Mario Brameshuber	2015_10_14  v1.2 updated  [Matlab 7.1 (R14)]
%   contact: mbrameshuber@yahoo.com
%   check fot program updates @ http://biophysics.iap.tuwien.ac.at/
%==========================================================================

a = MAT';
ind = size(a,2);
vec = 1:size(a,2);

for i = 1:ind*fraction
    r = ceil(ind*rand);
    vec(r) = [];
    ind = length(vec);
end


MAT = a(vec);
MAT = MAT';