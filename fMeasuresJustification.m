function fMeasuresJustification
% 2018 June 13
% Qi Wang


%% similar as fPrepareMeasures.mat)
rawdata = readtable('measures2015.csv', 'ReadVariableNames',true);
Dat0 = rawdata{:,2:end};
Dat = readtable('measures_current2015.csv', 'ReadVariableNames',true);
Dat = Dat{:,2:end};

measureList_new = readtable('list_2018_06_05.csv');
[subdomains_new, ia, ic] = unique(measureList_new.subdomain,'stable');

m0 = length(subdomains_new); % num subdomain
for j0=1:m0                                       % for each subdomain
    I_DS0{j0} = find(ic==j0);                    % get measures
    n_DS0(j0) = length(I_DS0{j0});                % number
end
% I_DS0 % coresponding measures col

n0 = length(measureList_new.subdomain);
for i0=1:n0                                       % for each current measure
    DS0(i0,:)  = string([ measureList_new{i0,3}{1}(1:2) ' ' measureList_new{i0,4}{1}(1:3) ]); % domain/subdomain abbreviation
    ttl0{i0}   = measureList_new{i0,2}{1};                   % measure title
end
[ u_DS0,~,iu0 ] = unique(DS0,'stable');             % 17 domains/subdomains

direc0 = strcmp(measureList_new.direction, '-');% positive or negative
direc0 = direc0';
Direc0 = zeros(size(direc0));               % change or not

%% Choose optimal measures
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
% keyboard

%% Factor loading
%----------------------------------------------------------------------------------------------
% CFA
%----------------------------------------------------------------------------------------------
fprintf('\nCFA of original measures\n')
for j0=1:m0                                           % for each current subdomain
    if n_DS0(j0)>1                                    % if at least two measures
        Is0   = I_DS0{j0};                            % collect the measures in DS0
        CX    = nancov(Dat0(:,Is0));                  % get the covariance
        lam   = fCFA( CX );                           % get the loadings
        %         lam   = lam/(0.5/3);                          % scale to 1
        lam   = lam/(sum(abs(lam)));                  % scale to 1
        CA    = fCAlpha(Dat0(:,Is0),Direc0max(Is0));  % default Cronbach's alpha
        fprintf('\n%s (%.2f)\n',u_DS0(j0,:),CA)
        %         lam =  lam*sign(lam(1));
        [~, maxlam_indx] = max(abs(lam));
        lam =  lam*sign(lam(maxlam_indx));
        for I0=1:n_DS0(j0)
            i0 = Is0(I0);
            %             if strcmp(Signs0(i0),'-'), lam(I0) = -lam(I0); end
            fprintf('  %s (%s) %5.2f %s\n',Signs0(i0),Signs(i0),lam(I0),ttl0{i0})
        end
        LAM{j0} = lam;
        
    end
end

% keyboard

%% Choose core subsets
for j0=1:m0                                           % for each current subdomain
    if n_DS0(j0)>1                                    % if at least two measures
        lam = LAM{j0};
        core_measures_index{j0} = [];
        Is0   = I_DS0{j0};
        for I0=1:n_DS0(j0)
            i0 = Is0(I0);
            if strcmp(Signs0(i0), Signs(i0)) & abs(lam(I0)) > 0.15
                core_measures_index{j0} = [core_measures_index{j0} Is0(I0)];
            end
        end
        
    end
end



% keyboard

%% Combine with the current measures
measureList_c = readtable('list_2016_12_13.csv');
direct0_current = strcmp(measureList_c.direction, '-');
for j0=1:m0                                           % for each current subdomain
    non_core_index{j0} = setdiff(I_DS0{j0},core_measures_index{j0});
end


for j0=1:m0                                           % for each current subdomain
    temp = strcmp(measureList_c.subdomain, subdomains_new{j0});
    cur_index{j0} = temp;
    %     datsub_cur = Dat(:,temp);
    datsub_core = Dat0(:,core_measures_index{j0});
    %     datsub = [datsub_core datsub_cur];
    
    
    if core_measures_index{j0}
        minchoose = 1;
        count = 1;
        index = find(temp);
        data_noncore = Dat(:,non_core_index{j0});
        data_current = Dat0(:,index);
        data_all = [data_noncore data_current];
        signs_noncore = direc0(:,non_core_index{j0});
        signs_current = direct0_current(index);
        signs_all = [signs_noncore'; signs_current];
        [~, numall] = size(data_all);
        numchoose = (1:numall)';
        CAmax_temp = [];
        Signs_temp = [];
        includeSavedtemp = zeros(numall,1);
        clear ID mea_domain covmatrix corre Id_signs
        [CAmax_temp(count), Signs_temp{count}] = fCAmaxperm_core(datsub_core, []);
    else
        minchoose = 2;
        count = 0;
        index = find(temp);
        data_noncore = Dat(:,non_core_index{j0});
        data_current = Dat0(:,index);
        data_all = [data_noncore data_current];
        signs_noncore = direc0(:,non_core_index{j0});
        signs_current = direct0_current(index);
        signs_all = [signs_noncore'; signs_current];
        [~, numall] = size(data_all);
        numchoose = (1:numall)';
        num_sub_cur = length(index);
        CAmax_temp = [];
        Signs_temp = [];
        includeSavedtemp = [];
        clear ID mea_domain covmatrix corre Id_signs
    end
    
    
    
    for kk = minchoose:min(4,length(numchoose))
        combinations = nchoosek(numchoose, kk);
        %         MA = zeros(length(numchoose),length(combinations));
        %         MA(combinations') = 1;
        
        [nrow, ~] = size(combinations);
        for cc = 1:nrow
            count = count + 1;
            datsub_cur = data_all(:, combinations(nrow, :));
            %             datasub_noncore = Dat0(:,)
            datsub = [datsub_core datsub_cur];
            clear ID mea_domain covmatrix corre Id_signs
            [ CAmax_temp(count), Signs_temp{count}] = fCAmaxperm_core(datsub_core, datsub_cur);
            includeSavedtemp = [includeSavedtemp ismember(numchoose, combinations(cc, :))];
        end
    end
    [nr, nc] = size( includeSavedtemp);
    [sortCA, sortCA_idx] = sort(CAmax_temp,'descend');
    CAsSaved{j0} = sortCA(1:min(nc, 20));
    SIGNS_optimal{j0} = Signs_temp;
    includeSavedtemp = includeSavedtemp(:, sortCA_idx(1:min(nc, 20)));
    %     includeSavedtemp = [ones(length(index),min(nc, 20));includeSavedtemp];
    includeSaved{j0} = includeSavedtemp;
    
    Signs_temp = Signs_temp(:, sortCA_idx(1:min(nc, 20)));
    SignsSaved{j0} = Signs_temp;
    Signs_REAL = [];
    for j = 1:min(nc, 20)
        inclu_temp = includeSavedtemp(:,j);
        signs_change = Signs_temp{j};
        if isempty(signs_change)
            Sign_real = zeros(size(inclu_temp));
        else
            Sign_real = zeros(size(inclu_temp));
            loc = find(inclu_temp ~=0);
            for jj = 1:length(loc)
                IN = loc(jj);
 
                if signs_change(jj)
                    Sign_real(IN) = 1-signs_all(IN);
                else
                    Sign_real(IN) = signs_all(IN);
                end
            end
            
        end
        Signs_REAL(:,j) = Sign_real;
    end
    SignsSaved{j0} = Signs_REAL;
    
    
    
    
    
    %
    %
    %     if (sum(Direc0max(Is0)) > length(direc0(Is0)) / 2)
    %         Direc0max(Is0) = 1 - Direc0max(Is0);
    %     end
    %     for i0=Is0'
    %         %                 fprintf('  %s (%s) %s\n',signs(direc0max(i0)+1),signs(direc0(i0)+1),ttl0{i0})
    %
    %         direc0max(i0) = mod(Direc0max(i0)+direc0(i0),2);  % used only for output
    %
    %         %                 fprintf('  %s (%s) %s\n',signs(direc0(i0)+1),signs(direc0max(i0)+1),ttl0{i0})
    %     end
    %
    %     CAlast = CAmax;
    %     Signs  = signs(direc0max+1);
    %     Signs0 = signs(direc0+1);
    %     if flag.Direc0, Direc0max = Direc0;
    %         CAlast = CA;
    %     end
    %     signs_all = direc0max(Is0);
    %
    %     mea_domain = Dat0(:,Is0).*(-1).^Direc0max(Is0)+ones(length(Dat0(:,Is0)),1)*(Direc0max(Is0)) ;
    %     dataplot {j0} = mea_domain;
    %     CA_plot{j0} = CAlast;
    %     [~,num] = size(mea_domain);
    %     %         for nid = 1: num
    %     %             Id_signs{nid} = [signs(signs_all(nid)+1) num2str(current_ID(Is0(nid)))];          % add the signs in the name of measure id
    %     %
    %     %         end
    %     Letters = {'A','B','C','D','E','F','G','H','I','J','K'};
    %
    %     % covariance
    %     covmatrix = cov(mea_domain,'omitrows');
    %     % correlation
    %     corre = corrcov(covmatrix);
    %
    
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

function [ CAmax,direcmax ] = fCAmaxperm_core(core, peri)
%----------------------------------------------------------------------------------------------
% max Cronbach's alpha
%----------------------------------------------------------------------------------------------
[n_core, ~] = size(core);
if ~isempty(peri)
    n_dat  = size(peri,2);                          % number of columns
    n_perm = 2^(n_dat-1);                          % number of permutations
    CAmaxs = zeros(n_dat,1);                       % initialize
    direcs = zeros(n_perm,n_dat);
    for jperm=1:n_perm                             % for all permutations
        direcs(jperm,:) = decimalToBinaryVector(jperm-1,n_dat);
        CAmaxs(jperm)   = fCAlpha([core peri],[zeros(1, n_core) direcs(jperm,:)]);
    end
    [ CAmax,Jmax ] = max(CAmaxs);                  % get the maximum
    direcmax       = direcs(Jmax,:);               % get the signs
else
    CAmax = fCAlpha(core, zeros(1, n_core));
    direcmax = [];
end


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


