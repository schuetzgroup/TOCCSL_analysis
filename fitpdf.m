%**************************************************************************
% Estimates the fraction of n-mers for a given monomer and multimer
% brightness distribution
%**************************************************************************
function OUT1 = fitpdf(single,multi,lim,fakt,noc,showfig);
%==========================================================================
% Call:   OUT1 = fitpdf(single,multi,lim,fakt,noc,showfig)
%         example:  fitpdf(monomers, multimers, 2000,2,4)
%
% Input:  single:  integrated intensity vector of MONOMERS (single)
%         multi:   integrated intensity vector of MULTIMERS (multi)
%                  recommended unit: detected photons (not counts)
%         lim:     intensity limit (time consuming convolution!)
%         fakt:    correction factor for error of brightness estimates /
%                  smoothing factor (choose e.g. 2)
%         noc:     number of convolutions (expected "noc"-mers)
%         showfig: 0 for no figure display, 1 or no argument for displaying 
%                  figure
%
% Output: probability density function, pdf 
%         Result:  fraction of monomers, dimers, etc.
%         ci: confidence intervals
%	
%   changes: v1.2   sum of all fitted x-mer fractions is one
%
%   Mario Brameshuber	2015_10_14  v1.2 updated  [Matlab 7.1 (R14)]
%   contact: mbrameshuber@yahoo.com
%   check fot program updates @ http://biophysics.iap.tuwien.ac.at/
%==========================================================================

%--------------------------------------------------------------------------
%Initialization
if size(single,2) > 5; single = single(:,5);end
if size(multi,2) > 5; multi = multi(:,5);end
if nargin < 6, showfig = 1; end
x = 0:1:lim;
y = zeros(1,max(size(x)));
y1 = y;

%--------------------------------------------------------------------------
%Calculate monomer pdf

for i = 1:max(size(single))
    mju = single(i);
    sig = fakt*sqrt(single(i));
    y1 = y1+normpdf(x,mju,sig);
end

y1=y1/(x(2)*sum(y1));   %normalization

if showfig == 1
disp(' '); disp('Monomer data:');
%maximum value:
disp(['max value:    ',num2str(x(find(y1==max(y1)))),' counts']);
%mean value:   mean=integrate(x*pdf)
disp(['mean value:   ',...
    num2str(sum(y1.*x*(max(x)-min(x))/max(size(x))),4),' counts']);
%median value:
disp(['median value: ',num2str(median(single),4),' counts']);
end
%standard deviation:    std=sqrt(integrate(x^2*pdf)-(integrate(x*pdf))^2)
%                       std=sqrt(integrate(x^2*pdf)-me^2)
%std_single=sqrt(sum(y1.*x.^2*(max(x)-min(x))/max(size(x)))-mean_single^2)

%--------------------------------------------------------------------------
%Calculate pdf of overall distribution

for i=1:max(size(multi))
    mju=multi(i);
    sig=fakt*sqrt(multi(i));
    y=y+normpdf(x,mju,sig);
end

y=y/(x(2)*sum(y));  %normalization

if showfig == 1
disp(' '); disp('Cluster data:');
%maximum value:
disp(['max value:    ',num2str(x(find(y==max(y)))),' counts']);
%mean value:   mean=integrate(x*pdf)
disp(['mean value:   ',...
    num2str(sum(y.*x*(max(x)-min(x))/max(size(x))),4),' counts']);
%median value:
disp(['median value: ',num2str(median(multi),4),' counts']);
end
%standard deviation:    std=sqrt(integrate(x^2*pdf)-(integrate(x*pdf))^2)
%                       std=sqrt(integrate(x^2*pdf)-me^2)
%std_multi=sqrt(sum(y.*x.^2*(max(x)-min(x))/max(size(x)))-mean_multi^2)

%--------------------------------------------------------------------------
%Calculation of n-mer distribution by convolution integrals

%set number of convolutions
if nargin < 5;noc=2;end

%calculate ditributions for 1:noc colocalized fluorophores
zz=convolution(y1,noc);
z=[];
for i=1:noc
    z=[z;zz{i}/sum(sum(zz{i}))];
end
y=[y,zeros(1,max(size(z))-max(size(y)))];
x=0:1:(max(size(y))-1);
y1=z(1,:);

%##########################################################################
%initialization
options = optimset('lsqcurvefit');
options.MaxFunEvals=options.MaxFunEvals*1000;
options.Display='off';
para=zeros(noc,1);
fun1=inline('para*z','para','z');
%%%% fun1=inline('para*z(1:end-1,:)+z(end,:)*(1-sum(para))','para','z');
%%%% sfit=0.5.^(1:1:noc-1);
sfit=0.5.^(1:1:noc);

%1.) fitting with nlinfit
%Fast but allows negative paramters
%[OUT1,resid,J]=nlinfit(z,y,fun1,sfit);

%2.) lsqcurvefit:
%%%% lowerb=zeros(noc-1,1); upperb=ones(noc-1,1);
lowerb=zeros(noc,1); upperb=ones(noc,1);
[OUT1,rn,resid,ef,op,lm,J]=lsqcurvefit(fun1,sfit,z,y,lowerb,upperb,options);
J=full(J); %conversion of a sparse to a full matrix
%OUT1(2)=1-OUT1(1);

%3.) lsqnonneg: z*OUT=y
%Fast but calculates no Jacobian
%OUT2=lsqnonneg(z',y')
if showfig == 1
disp(' '); disp(['Result: ',num2str(OUT1,3)]);
%calculate confidence intercals for fitting paramters:
ci=nlparci(OUT1,resid,J),ferr=(ci(:,2)-ci(:,1))';
disp(['-> fit errors: ',num2str((ci(:,2)-ci(:,1))',2)]);
end

%##########################################################################

if showfig>0
start=ones(noc,1);

figure; hold on; plot(x,y,'k','LineWidth',1.5)

if max(size(start))==20
    for i=1:15
        if i<=8
            plot(x,OUT1(i)*z(i,:),'b','LineWidth',1.5)
        else
            plot(x,OUT1(i)*z(i,:),'b','LineWidth',1)
        end
    end
else
    for i=1:noc
        if i<=6
            plot(x,OUT1(i)*z(i,:),'b','LineWidth',1.5)
        else
            plot(x,OUT1(i)*z(i,:),'b','LineWidth',1)
        end
    end
end

%plot fitted pdf
zneu=zeros(size(z(1,:)));
for i=1:noc
    zneu=zneu+z(i,:)*OUT1(i);
end
plot(x,zneu,'r')

%title('fitpdf','FontSize',16);
xlabel('brightness (counts)','FontSize',16);
ylabel('pdf (counts^-^1)','FontSize',16);
set(gca,'Box','off','LineWidth',1.5,'FontSize',14,'XLim',[0,lim]);
%% Create legend
switch noc
    case 2
        legend({['cluster signal (',num2str(size(multi,1)),'dp)'],...
            ['monomers (',num2str(size(single,1)),'dp)'],...
            ['dimers (',num2str(round(OUT1(2)*100)),'%)'],...
            'fit'},'FontSize',14);
    case 3
        legend({['cluster signal (',num2str(size(multi,1)),'dp)'],...
            ['monomers (',num2str(size(single,1)),'dp)']...
            ,['dimers (',num2str(round(OUT1(2)*100)),'%)'],...
            ['trimers (',num2str(round(OUT1(3)*100)),'%)'],'fit'},...
            'FontSize',14);
    case 4
        legend({['cluster signal (',num2str(size(multi,1)),'dp)'],...
            ['monomers (',num2str(size(single,1)),'dp)'],...
            ['dimers (',num2str(round(OUT1(2)*100)),'%)'],...
            ['trimers (',num2str(round(OUT1(3)*100)),'%)'],...
            ['tetramers (',num2str(round(OUT1(4)*100)),'%)'],'fit'},...
            'FontSize',14);
end
hold off
end
%--------------------------------------------------------------------------
