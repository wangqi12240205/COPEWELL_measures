function fMeasuresJustification
% 2018 June 13
% Qi Wang


%% similar as fPrepareMeasures.mat)
rawdata = readtable('measures2015.csv', 'ReadVariableNames',true);
Dat0 = rawdata{:,2:end};
Dat = readtable('measures_current2015.csv', 'ReadVariableNames',true);
Dat = Dat{:,2:end};

measureList = readtable('list_2018_06_05.txt');
[subdomains, ia, ic] = unique(measureList.subdomain,'stable');
m0 = length(subdomains); % num subdomain
for j0=1:m0                                       % for each subdomain
    I_DS0{j0} = find(ic==j0);                    % get measures
end
% I_DS0 % coresponding measures col

n0 = length(measureList.subdomain);
for i0=1:n0                                       % for each current measure
    DS0(i0,:)  = char([ measureList{i0,3}{1}(1:2) ' ' measureList{i0,4}{1}(1:3) ]); % domain/subdomain abbreviation
    ttl0{i0}   = measureList{i0,2}{1};                   % measure title
end
[ u_DS0,~,iu0 ] = unique(DS0,'rows','stable');             % 17 domains/subdomains

direc0 = strcmp(measureList.direction, '-');% positive or negative
direc0 = direc0';
Direc0 = zeros(size(direc0));               % change or not

flag.maxdirec = 1;

for tempDiret= 0:1
    flag.Direc0 = tempDiret;
    if flag.maxdirec
        %----------------------------------------------------------------------------------------------
        % Summary of measures
        %----------------------------------------------------------------------------------------------
        signs               = ['+','-'];
        Direc0max = Direc0;
        for j0=1:m0                                       % for each current subdomain
            clear ID mea_domain covmatrix corre Id_signs
            Is0   = I_DS0{j0};                            % collect the measures in DS0
            CA    = fCAlpha(   Dat0(:,Is0),Direc0(Is0));  % default Cronbach's alpha
            [ CAmax,Direc0max(Is0) ] = fCAmaxperm(Dat0(:,Is0));
            %   fprintf('%s (%i) %5.2f (%5.2f)\n',u_DS0(j0,:),n_DS0(j0),CA,CAmax)
            if (sum(Direc0max(Is0)) > length(direc0(Is0)) / 2)
                Direc0max(Is0) = 1 - Direc0max(Is0);
            end
            for i0=Is0'
                %                 fprintf('  %s (%s) %s\n',signs(direc0max(i0)+1),signs(direc0(i0)+1),ttl0{i0})
                
                direc0max(i0) = mod(Direc0max(i0)+direc0(i0),2);  % used only for output
                
                %                 fprintf('  %s (%s) %s\n',signs(direc0(i0)+1),signs(direc0max(i0)+1),ttl0{i0})
            end
            
            CAlast = CAmax;
            Signs  = signs(direc0max+1);
            Signs0 = signs(direc0+1);
            if flag.Direc0, Direc0max = Direc0;
                CAlast = CA;
            end
            signs_all = direc0max(Is0);
            
            mea_domain = Dat0(:,Is0).*(-1).^Direc0max(Is0)+ones(length(Dat0(:,Is0)),1)*(Direc0max(Is0)) ;
            dataplot {j0} = mea_domain;
            CA_plot{j0} = CAlast;
            [~,num] = size(mea_domain);
            %         for nid = 1: num
            %             Id_signs{nid} = [signs(signs_all(nid)+1) num2str(current_ID(Is0(nid)))];          % add the signs in the name of measure id
            %
            %         end
            Letters = {'A','B','C','D','E','F','G','H','I','J','K'};
            
            % covariance
            covmatrix = cov(mea_domain,'omitrows');
            % correlation
            corre = corrcov(covmatrix);
            
        end
    end
    if flag.Direc0    % original
        save('dataplot_original.mat','dataplot','CA_plot','u_DS0');
    else
        save('dataplot_currentmaxv2.mat','dataplot','CA_plot','u_DS0');
    end
end
keyboard







function [ lambda,r2jj ] = fCFA( CX, niter )
%----------------------------------------------------------------------------------------------
% CFA
%----------------------------------------------------------------------------------------------
CXjj = diag(CX);                                 % cX,jj
CXh  = diag(sqrt(1./CXjj));                      % cX,jj^-1/2
cX   = CXh*CX*CXh;                               % rhoX
cXi  = inv(cX);                                  % cX^-1
r2jj = ones(size(CXjj));                         % initialize rho,jj^2
if nargin<2, niter = 100; end
for iter=1:niter
    ceps     = diag(r2jj);                       % cepsilon
    [lam,om] = eigs(cXi,ceps,1,'SM');            % smallest eigenvector and eigenvalue
    scale    = sqrt((1 - om)/(lam'*cXi*lam));
    lambda   = scale*lam;                        % scaled lambda
    r2jj     = 1 - lambda.^2;
end
r2jj   = r2jj.*CXjj;
lambda = CXh\lambda;


function [ CAmax,Imax,sgn ] = fCAmax( CA, dat0, dat, I_match, Direc0 )
%----------------------------------------------------------------------------------------------
% Find maximum CA
%----------------------------------------------------------------------------------------------
CAs   = zeros(2,size(dat,2));
for i=I_match                                    % for each original measure
    for k=1:2                                    % +/- sign of the measure
        CAs(k,i) = fCAlpha( [dat0 dat(:,i)],[Direc0 k-1]);
    end
end
if size(dat0,2)==1 && CA==1 && ~isempty(I_match) % if only 1 measure (and it is 1)
    CA0 = 0;                                     % always add original measure
else
    CA0 = CA;
end
[ CAmax,IKmax ] = max([CAs(:)]);           % take the maximum alpha
Imax = ceil(IKmax/2);                            % get the measure index (n+1 if CA is max)
Kmax = IKmax - 2*Imax + 2;                       % +/- sign for Kmax = 1,2
if Kmax==2, sgn = 1; else sgn = 0; end

function [CAmaxComb, combIndex] = fCAmaxCombination(Dat,dirction)
[~,numMea] = size(Dat);
CAmaxComb = 0;
combIndex =0;
if (numMea < 2)
    CAmaxComb = 1;
    combIndex = 1;
else
    for ii = 2:numMea
        C = combnk(1:numMea,ii);
        [numComb, ~] = size(C);
        for jj = 1:numComb
            datComb = Dat;
            datComb(:,C(jj,:)) = 2*Dat(:,C(jj,:));
            CA = fCAlpha(datComb, dirction);
            if (CA>CAmaxComb)
                CAmaxComb = CA;
                combIndex = C(jj,:);
            end
            
        end
    end
end


function [ CAmax,direcmax ] = fCAmaxperm(dat)
%----------------------------------------------------------------------------------------------
% max Cronbach's alpha
%----------------------------------------------------------------------------------------------
n_dat  = size(dat,2);                          % number of columns
n_perm = 2^(n_dat-1);                          % number of permutations
CAmaxs = zeros(n_dat,1);                       % initialize
direcs = zeros(n_dat,n_dat);
for jperm=1:n_perm                             % for all permutations
    direcs(jperm,:) = decimalToBinaryVector(jperm-1,n_dat);
    CAmaxs(jperm)   = fCAlpha(dat,direcs(jperm,:));
end
[ CAmax,Jmax ] = max(CAmaxs);                  % get the maximum
direcmax       = direcs(Jmax,:);               % get the signs


function bin = decimalToBinaryVector(dec,num_digits)
%----------------------------------------------------------------------------------------------
% decimal to binary
%----------------------------------------------------------------------------------------------
Bin = dec2bin(dec);                               % convert to binary string
M   = length(Bin);                                % number of binary digits from dec
if nargin==1
    num_digits = M;                               % if num_digits not specified, use M
end
bin = zeros(1,num_digits);                        % initialize bin

for i=1:num_digits
    if i<=M && Bin(M-i+1) == '1'                  % if within M
        bin(i) = 1;                               % set
    end
end


function alp = fCAlpha(Y,Direc0)
%----------------------------------------------------------------------------------------------
% Cronbach's alpha
%----------------------------------------------------------------------------------------------
global CAcount

CAcount = CAcount + 1;

[ ~,K ] = size(Y);                                 % K = number of components
Ikeep   = ~isnan(sum(Y,2));                        % get good rows
Y       = Y(Ikeep,:);                              % keep only good rows
if K==1, alp = 1; return, end                      % default
for k=1:K
    if Direc0(k), Y(:,k) = 1-Y(:,k); end           % change sign, as needed
    sig(k) = nanvar(Y(:,k));                       % variance of each score
end
X       = nanmean(Y,2)*K;                          % total scores
alp     = K / (K-1) * ( 1 - sum(sig)/nanvar(X) );  % Cronbach's alpha
if alp>1, keyboard, return, end

