%**************************************************************************
% calculate convolution integrals
%**************************************************************************
function y = convolution(z,d)
%==========================================================================
% Call:   y = convolution(z,d)
%
% Input:  monomer pdf (z), number of convolutions (d)
%        
% Output: multimer pdfs in cell array; y{1} = monomer pdf, y{2} = dimer pdf
%    
%        
%   Mario Brameshuber	2015_10_14  v1.2 updated  [Matlab 7.1 (R14)]
%   contact: mbrameshuber@yahoo.com
%   check fot program updates @ http://biophysics.iap.tuwien.ac.at/
%==========================================================================

%calculate individual pdfs
y{1} = z;
for i = 2:d
    z = conv(y{1},y{i-1});
    y{i} = z/sum(z);
end
    
%normalize pdf lengths
z2 = max(size(y{d}));

for i = 0:d-2
    z1 = y{d-i-1};
    z1 = [z1,zeros(1,z2-max(size(z1)))];
    y{d-i-1} = z1;
end