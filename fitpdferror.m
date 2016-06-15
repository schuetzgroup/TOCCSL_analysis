%**************************************************************************
% Error estimation for fitpdf by drawing random 50% sub-samples out of
% monomer and multimer distribution
%**************************************************************************
function result = fitpdferror(single,multi,lim,fakt,noc,runs);
%==========================================================================
% Call:   result = fitpdferror(single,multi,lim,fakt,noc,runs)
%         example:  fitpdferror(monomers, multimers, 2000,2,4,100)
%
%         for detailed description refer to program fitpdf
%
% Input:  runs: number of iterations
% Output: n-mer fractions for every run
%         mean/std/S.E.M. of n-mer fractions
%	
%   Mario Brameshuber	2015_10_14  v1.2 updated  [Matlab 7.1 (R14)]
%   contact: mbrameshuber@yahoo.com
%   check fot program updates @ http://biophysics.iap.tuwien.ac.at/
%==========================================================================

%--------------------------------------------------------------------------
%Initialization
if size(single,2) > 5; single = single(:,5);end
if size(multi,2) > 5; multi = multi(:,5);end

h = waitbar(0,'Please wait...');

result = [];

for i = 1:runs
    
    out = fitpdf(subset(single,0.5),subset(multi,0.5),lim,fakt,noc,0);
    %out = fitpdf(single,subset(multi,0.5),lim,fakt,noc,0);
    result = [result;out];
    waitbar(i / runs);
    
end

result.res = result;
result.mean_std_SEM = [mean(result.res); std(result.res); std(result.res)/sqrt(2)];
close(h) 