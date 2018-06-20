function fOptimizeMeasuresv2
%----------------------------------------------------------------------------------------------
% Prepare the data
%----------------------------------------------------------------------------------------------
global flag n n0 cor cor0 cor1 DS0 DS ttl0 ttl Dat0 Dat Fips current_ID
global m m0 u_DS u_DS0 u_Dom u_Dom0 u_Sub u_Sub0 n_DS n_DS0 direc0 Direc0 I_DS I_DS0 alp0
global CAcount

flag.load           = 1;                          % 1-load data (default), 0-read data

fReadMeasures                                     % load or read data

flag.summary        = 1;
flag.crosscor_print = 0;
flag.cor_plot       = 0;
flag.measures_corr  = 0;

flag.append_one     = 0;
flag.subtract_one   = 0;

CAcount = 0;

if flag.summary
    %----------------------------------------------------------------------------------------------
    % Summary of measures
    %----------------------------------------------------------------------------------------------
    flag_print  = 0;
    I_DS0_match = fDS_match( flag_print );
    flag.p              = 1;                          % print signs
    flag.Direc0         = 0;                          % use original direc
    signs               = ['+','-'];
    Direc0max = Direc0;
    for j0=1:m0                                       % for each current subdomain
        Is0   = I_DS0{j0};                            % collect the measures in DS0
        CA    = fCAlpha(   Dat0(:,Is0),Direc0(Is0));  % default Cronbach's alpha
        [ CAmax,Direc0max(Is0) ] = fCAmaxperm(Dat0(:,Is0));
        %   fprintf('%s (%i) %5.2f (%5.2f)\n',u_DS0(j0,:),n_DS0(j0),CA,CAmax)
        if flag.p
            if (sum((mod(Direc0max(Is0)+direc0(Is0),2))) > length(direc0(Is0)) / 2)
                Direc0max(Is0) = 1 - Direc0max(Is0);
            end
            for i0=Is0'
                %                 fprintf('  %s (%s) %s\n',signs(direc0max(i0)+1),signs(direc0(i0)+1),ttl0{i0})
                
                direc0max(i0) = mod(Direc0max(i0)+direc0(i0),2);  % used only for output
                
                %                 fprintf('  %s (%s) %s\n',signs(direc0(i0)+1),signs(direc0max(i0)+1),ttl0{i0})
            end
        end
    end
    Signs  = signs(direc0max+1);
    Signs0 = signs(direc0+1);
    if flag.Direc0, Direc0max = Direc0; end
    
    %----------------------------------------------------------------------------------------------
    % Replace one measure
    %----------------------------------------------------------------------------------------------
    % fprintf('\nReplace one current measure with original measure\n')
    fprintf('\nCURRENT MEASURES: Substracting some negative current measures and adding some original measures\n')
    for j0=1:m0                                              % for each current subdomain
        clear ID mea_domain covmatrix corre Id_signs I_DS0_match_temp
        if n_DS0(j0)>2                                       % if at least two measures
            Is0   = I_DS0{j0};                               % collect the measures in DS0
            CA    = fCAlpha(Dat0(:,Is0),Direc0max(Is0));     % default Cronbach's alpha
            CAmaxs= zeros(n_DS0(j0),1); Imaxs = CAmaxs; sgns = CAmaxs;
            Is0sub= Is0;
            I_DS0_match_temp = I_DS0_match{j0};
            I_DS0_add  =[];
            sgns_add   = [];                                 % initialiaze                          
            fprintf('\nSubdomain: %s\nCA: %5.2f',u_DS0(j0,:),CA)
            CAlast = fCAlpha([Dat0(:,Is0sub) Dat(:,I_DS0_add)],[Direc0max(Is0sub) sgns_add]);
            %----------substract current measure(-)
            for j0sub=1:n_DS0(j0)                            % subtract current measure
                if (Signs(Is0(j0sub)) == '-')
                    temp = setdiff(Is0sub,Is0(j0sub));       % get reduced set
                    CAsubstract = fCAlpha([Dat0(:,temp) Dat(:,I_DS0_add)],[Direc0max(temp) sgns_add]);
                    if (CAsubstract < CAlast)                % if after substraction, CA decrease, then don't substract
                        CAsubstract = CAlast;
                        Is0sub = Is0sub;
                    else
                        Is0sub = temp;                       % if original set is empty
                        fprintf('\nSubtract current measure ID %d: %5.2f',current_ID(Is0(j0sub)),CAsubstract);
                        if (length(I_DS0_match_temp) == 0)
                        fprintf('\nNothing to be added');
                        else
             %---------replace with the best original measure           
                        [ CAmaxs(j0sub),Imaxs(j0sub),sgns(j0sub) ] = ...
                            fCAmax( CAsubstract, [Dat0(:,Is0sub) Dat(:,I_DS0_add)], Dat, I_DS0_match_temp, [Direc0max(Is0sub) sgns_add] );
                        I_DS0_match_temp = setdiff(I_DS0_match_temp, Imaxs(j0sub));
                        I_DS0_add = [I_DS0_add; Imaxs(j0sub)];
                        sgns_add  = [sgns_add sgns(j0sub)];
                        CAlast    = CAmaxs(j0sub);
                        fprintf('\nAdd original measure No.%d: %5.2f',Imaxs(j0sub),CAmaxs(j0sub));
                        end
                    end
                    
                    
                end
            end
            %---------Add original measure
            while (length(I_DS0_match_temp) ~= 0 && CAlast < 0.6)
                [ CAmaxs,Imaxs,sgns ] = ...
                    fCAmax( CAlast, [Dat0(:,Is0sub) Dat(:,I_DS0_add)], Dat, I_DS0_match_temp, [Direc0max(Is0sub) sgns_add] );
                I_DS0_match_temp = setdiff(I_DS0_match_temp, Imaxs);
                I_DS0_add = [I_DS0_add; Imaxs];
                sgns_add  = [sgns_add sgns];
                CAlast    = CAmaxs;
                fprintf('\nAdd original measure No.%d: %5.2f',Imaxs,CAmaxs);
            end
           
            
            Dat_comb =  [Dat0(:,Is0sub) Dat(:,I_DS0_add)];
            signs_all = [Direc0max(Is0sub) sgns_add];
            
            
            
            
            %---------scale measures by 2 
            ID = [Is0sub; I_DS0_add];
            [~,num] = size(Dat_comb);
            for nid = 1: num
                if (nid <= length(Is0sub))
                    Id_signs{nid} = [signs(signs_all(nid)+1) num2str(current_ID(ID(nid)))];          % add the signs in the name of measure id
                else
                    Id_signs{nid} = [signs(signs_all(nid)+1) 'Ori' num2str(ID(nid))]; 
                end
            end
            [CAmaxComb, combIndex] = fCAmaxCombination(Dat_comb,signs_all);
            num_combinde = length(combIndex);
            for i = 1: num_combinde
                id= Id_signs{combIndex(i)};
                fprintf('\nScale measure %s by 2',id)
            end
             fprintf('\nCA: %5.2f',CAmaxComb);
             fprintf('\nEnd\n');
             CAlast = CAmaxComb;
             mea_domain = Dat_comb;
             mea_domain(:,combIndex) = 2*Dat_comb(:,combIndex);
            
            % covariance
            covmatrix = cov(mea_domain,'omitrows');
            % correlation
            corre = corrcov(covmatrix);
            % plot the figure
            figure(500+j0),
            fig=gcf; fig.Position = [100 300 1300 500];                           % figure position and size
            subplot(1,2,1)                                                        % subplot (1,2,1)
            fig = gcf;
            fig.Name = u_DS0(j0,:);                                               % figure name is subdomain name
            
            clims = [-0.5:0.25:1];                                                % color label
            imagesc(corre,[-0.5 1]),colorbar('Ticks',clims)                       % image section with color bar
            title(['Directions maximizing alpha, \alpha = ' num2str(round(CAlast,3))],'Fontsize',18) % title
            set(gca,'xtick',[1:num]),set(gca,'xticklabel',Id_signs,'Fontsize',15);% xtick
            set(gca,'ytick',[1:num]),set(gca,'yticklabel',Id_signs,'Fontsize',15);% ytick
            
            subplot(1,2,2)                                                        % subplot(1,2,2)
            [~, AX, ~, ~, ~] = plotmatrix(mea_domain);                            % ax of figure
            set(AX,'XLim',[-0.2,1.2]);                                            % set xlim
            set(AX,'YLim',[-0.2,1.2]);                                            % set ylim
            
            set(AX,'YTick',[0 0.5 1]);                                            % set xtick
            set(AX,'XTick',[0 0.5 1]);                                            % set ytick
            %set the text of correlation inside the plot matrix
            for j = 1: num
                for mm = 1: num
                    if j ~=mm
                        if corre(j,mm)>0
                            co = 'red';
                        else
                            co = 'black';
                        end
                        text(AX(j,mm),0,0,num2str(round(corre(j,mm),2)),'color',co)
                    end
                end
            end
            
            title('Correlation Matrix','Fontsize',18)
            saveas(gca, ['figure_corr/' u_DS0(j0,:) '.fig']);
           
        end
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



function I_DS0_match = fDS_match( flag_print )
%----------------------------------------------------------------------------------------------
% Match original and current measure sets
%----------------------------------------------------------------------------------------------
global m0 m DS n_DS0 n_DS u_DS0 u_DS I_DS0 I_DS ttl0 ttl flag n

%----------------------------------------------------------------------------------------------
% Print measures
%----------------------------------------------------------------------------------------------
if flag_print
    fprintf('\nCURRENT MEASURES\n\n')
    for j0=1:m0
        fprintf('%2i %s\n',j0,u_DS0(j0,:))
    end
    
    fprintf('\nCURRENT MEASURES\n')
    for j0=1:m0
        fprintf('\n%2i %s\n',j0,u_DS0(j0,:))
        for i=1:n_DS0(j0)
            fprintf('    %s\n',ttl0{I_DS0{j0}(i)})
        end
    end
end

flag_individual = 1;
%----------------------------------------------------------------------------------------------
% Organize measures
%----------------------------------------------------------------------------------------------
if flag_individual
    x           = NaN;
    I_DS_match = [
        2 11 x 16 x x 14 x x x x 11 x x x 6 6 ...
        6 14 7 14 2 17 7 x 14 x x x x x 10 x 14 ...
        17 2 1 ...
        17 17 17 17 17 2 16 x 16 ...
        14 6 6 ...
        x x 4 x 11 11 11 11 x x 6 6 ...
        6 x 6 6 6 x ...
        x x x ...
        x ...
        x 16 x 14 ...
        10 ...
        2 7 ...
        12 ...
        12 x x x ...
        11 x x ...
        10 2 4 x x ...
        13 ...
        x ...
        x x x x x x ...
        17 ];
    if flag_print
        fprintf('\nORIGINAL MEASURES\n')
        for j=1:m
            fprintf('\n%s (%i)\n',u_DS(j,:),n_DS(j))
            for i=1:n_DS(j)
                ii = I_DS{j}(i);
                fprintf('%3i %3i %s\n',I_DS_match(ii),ii,ttl{ii})
            end
        end
        
    end
else
    DS_match  = {
        [  5 ]                                % AO
        [  5 ]                                % AP
        [ 12 ]                                % AR
        [ 17 ]                                % AS
        [  3 ]                                % CE
        [  4 ]                                % CF
        [  6 ]                                % CH
        [  2 ]                                % CJ
        [    ]                                % CL
        [  1 11 ]                             % CN
        [    ]                                % CP (11th row, only included if flag.keep_all)
        [ 10 ]                                % CT
        [ 11 ]                                % GB
        [ 12 ]                                % GG
        [ 12 ]                                % GL
        [ 13 ]                                % PE
        [ 13 ]                                % PR
        [ 13 ]                                % PW
        [ 13 ]                                % RO
        [ 13 ]                                % SO
        };
    if ~flag.keep_all
        for j=11:m, DS_match{j} = DS_match{j+1}; end
        DS_match{end} = NaN;
    end
end

DS0_match = cell(m0,1);
% fprintf('\nORIGINAL MEASURES (August 2014)\n')
for j0=1:m0                                          % for each current measure
    DS0_match{j0} = [];
    if flag_print, fprintf('\n%2i %s\n',j0,u_DS0(j0,:)), end
    if ~flag_individual
        I_DS0_match{j0} = [];
        for j=1:m                                        % for each original measure
            if ismember(j0,DS_match{j})                  % if included
                DS0_match{j0}   = [   DS0_match{j0} j ]; % add
                I_DS0_match{j0} = [ I_DS0_match{j0} I_DS{j} ];
                if flag_print, fprintf('      %s\n',u_DS(j,:)), end
            end
        end
    else
        I_DS0_match{j0} = find(I_DS_match==j0);
    end
    if flag_print
        for i=I_DS0_match{j0}
            %           fprintf('  %3i %s %s\n',i,DS(i,:),ttl{i})
            fprintf('  %s %s\n',DS(i,:),ttl{i})
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
    if Direc0(k), Y(:,k) = -Y(:,k); end            % change sign, as needed
    sig(k) = nanvar(Y(:,k));                       % variance of each score
end
X       = nanmean(Y,2)*K;                          % total scores
alp     = K / (K-1) * ( 1 - sum(sig)/nanvar(X) );  % Cronbach's alpha
if alp>1, keyboard, return, end


function fReadMeasures
%----------------------------------------------------------------------------------------------
% Prepare the data
%----------------------------------------------------------------------------------------------
global flag n n0 cor cor0 cor1 DS0 DS ttl0 ttl Dat0 Dat Fips current_ID
global m m0 u_DS u_DS0 u_Dom u_Dom0 u_Sub u_Sub0 n_DS n_DS0 direc0 Direc0 I_DS I_DS0 alp0

if flag.load, load fPrepareMeasures, return, end  % load data
%----------------------------------------------------------------------------------------------
% Read the data
%----------------------------------------------------------------------------------------------
T    = readtable('new_table_2010.csv');           % get original measures
T0   = readtable('measures2010.csv');             % get current  measures
T0   = T0(2:end,:);                               % get rid of first row
Ind  = readtable('CopeWELL indices Aug 7.xls');   % get original domains/subdomains
Ind0 = readtable('list_2016_12_13.csv');          % get current  domains/subdomains

%----------------------------------------------------------------------------------------------
% Set basic parameter and variables
%----------------------------------------------------------------------------------------------
current_ID = Ind0.id;
n       = 152;                                    % number of measures
n0      = length(Ind0{:,1});                      % number of current measures
[ N,~ ] = size(T);                                % number of counties
fips    = T{:,1};                                 % 3142 FIPS codes
for I=1:N+1                                       % there is one more FIPS in current version
    fips0(I,1) = str2num(T0{I,1}{1});             % get the FIPS codes
end
[ Fips,ifips,ifips0 ] = intersect(fips,fips0);    % find the common FIPS
T       = T (:,2:n+1);                            % get rid of FIPS, keep only n=152 columns
T0      = T0(:,2:end);
T       = T (ifips, :);                           % get original data for common FIPS only
T0      = T0(ifips0,:);                           % get current  data for common FIPS only
N       = length(Fips);                           % number of common FIPS
for i=1:n
    fprintf('%3i %s\n',i,Ind{i,6}{1})
end

%----------------------------------------------------------------------------------------------
% Extract good original data
%----------------------------------------------------------------------------------------------
flag.keep_all = 0;                                % keep all data
flag.Pres     = 1;                                % use Presidential election data
if flag.keep_all
    Ikeep   = 1:n;                                % keep all
else
    Ikeep   = setdiff(1:n,...
        [23,25,26,28,32:45,49:57,67:69,72:81,95,97,118,120,121,125:128,141,143,152,153]);
    if flag.Pres
        Ikeep(find(Ikeep==58)) = 51;
    end
end
Ind     = Ind(Ikeep,:);                           % keep only good columns
T       = T(:,Ikeep);                             % take a subset of the data
n       = length(Ikeep);                          % new number of measures

%----------------------------------------------------------------------------------------------
% Get original headings
%----------------------------------------------------------------------------------------------
for i=1:n                                         % for each measure
    DomEx(i,1) = Ind{i,1}{1};                     % domain letters (A,C,G,P,R,S)
    SubEx(i,1) = Ind{i,3}{1};                     % subdomain letters (B,E,...,W)
    ttlEx{i}   = Ind{i,6}{1};                     % titles
    fprintf('\n%s\n%s\n',T.Properties.VariableNames{i},ttlEx{i})
end

for i0=1:n0                                       % for each current measure
    DS0(i0,:)  = char([ Ind0{i0,3}{1}(1:2) ' ' Ind0{i0,4}{1}(1:3) ]); % domain/subdomain abbreviation
    ttl0{i0}   = Ind0{i0,2}{1};                   % measure title
    direc0(i0) = strcmp(Ind0{i0,5},'-');
end
Direc0= zeros(size(direc0));                      % current data is already scaled

DSEx  = [DomEx SubEx];                            % combine domain and subdomain letters
DS0   = char(DS0);                                % make it a character array

u_Dom  = unique(DomEx);                           %  6
u_Dom0 = unique(DS0(:,1:2),'rows');               % 5  domains
u_Sub  = unique(SubEx);                           % 14
u_Sub0 = unique(DS0(:,4:5),'rows');               % 17 subdomains
[ u_DS, ~,iu  ] = unique(DSEx, 'rows');           % 19 unique combinations of domain/subdomain
[ u_DS0,~,iu0 ] = unique(DS0,'rows');             % 17 domains/subdomains
m     = length(u_DS);                             % number of domain/subdomain combinations
m0    = length(u_DS0);                            % 17

%----------------------------------------------------------------------------------------------
% Get data
%----------------------------------------------------------------------------------------------
Dat   = zeros(N,n);                               % original data
Dat0  = zeros(N,n0);
IDat  = 0;                                        % index
IDat0 = 0;

for j=1:m                                         % for each subdomain
    I_DSj   = find(iu==j);                        % get measures (Excel data)
    n_DS(j) = length(I_DSj);                      % number
    fprintf('\n%s (%i)\n',u_DS(j,:),n_DS(j))
    I_DS{j} = IDat+1:IDat+n_DS(j);                % save indeces
    for i=1:n_DS(j)                               % for each measure
        IDat        = IDat + 1;                   % increment index
        IEx         = I_DSj(i);                   % Excel file index
        ttl{IDat}   = ttlEx{IEx};                 % transfer Excel data using ordered index
        DS(IDat,:)  = DSEx(IEx,:);
        fprintf('  %3i %s\n',I_DS{j}(i),ttl{IDat})
        Dat(:,IDat) = T{:,I_DSj(i)};              % save measure data
    end
end

for j0=1:m0                                       % for each subdomain
    I_DS0{j0} = find(iu0==j0);                    % get measures
    n_DS0(j0) = length(I_DS0{j0});                % number
    fprintf('\n%s (%i)\n',u_DS0(j0,:),n_DS0(j0))
    for i0=1:n_DS0(j0)                            % for each measure
        fprintf('  %3i %s\n',I_DS0{j0}(i0),ttlEx{I_DS0{j0}(i0)})
        IDat0         = IDat0 + 1;                % increment index
        Dat0(:,IDat0) = T0{:,I_DS0{j0}(i0)};      % save measure data
    end
end

%----------------------------------------------------------------------------------------------
% Compute correlations, alpha
%----------------------------------------------------------------------------------------------
for i=1:n                                         % for each original measure
    for i0=1:n0                                   % for each current measure
        cor(i,i0)      = corr(Dat (:,i),Dat0(:,i0),'rows','pairwise'); % cross-correlation
        if i<=n0                                  % note that n0 < n
            cor0(i,i0) = corr(Dat0(:,i),Dat0(:,i0),'rows','pairwise'); % current correlation
        end
    end
    for j=1:n                                     % for each original measure
        cor1(i,j)      = corr(Dat (:,i),Dat (:,j), 'rows','pairwise'); % original correlation
    end
end

for j0=1:m0                                       % for each original subdomain
    alp0(j0) = fCAlpha(Dat0(:,I_DS0{j0}),Direc0(I_DS0{j0})); % Cronbach's alpha
end

save fPrepareMeasures
