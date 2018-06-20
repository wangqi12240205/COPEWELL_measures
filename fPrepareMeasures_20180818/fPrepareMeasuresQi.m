function fPrepareMeasuresQi
%----------------------------------------------------------------------------------------------
% Prepare the data
%----------------------------------------------------------------------------------------------
global flag n n0 cor cor0 cor1 DS0 DS ttl0 ttl Dat0 Dat Fips
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
    flag.p              = 1;                          % print signs
    flag.Direc0         = 0;                          % use original direc
    signs               = ['+','-'];
    fprintf('\nSummary of measures, Cronbach''s alpha\n')
    Direc0max = Direc0;
    fprintf('\nCURRENT MEASURES: Optimizing directionality (parentheses for Cronbach''s alpha and direction for revised directions)\n')
    for j0=1:m0                                       % for each current subdomain
        Is0   = I_DS0{j0};                            % collect the measures in DS0
        CA    = fCAlpha(   Dat0(:,Is0),Direc0(Is0));  % default Cronbach's alpha
        [ CAmax,Direc0max(Is0) ] = fCAmaxperm(Dat0(:,Is0));
        %   fprintf('%s (%i) %5.2f (%5.2f)\n',u_DS0(j0,:),n_DS0(j0),CA,CAmax)
        fprintf('\n%s %5.2f (%5.2f)\n',u_DS0(j0,:),CA,CAmax)
        if flag.p
            if (sum((mod(Direc0max(Is0)+direc0(Is0),2))) > length(direc0(Is0)) / 2)
                Direc0max(Is0) = 1 - Direc0max(Is0);
            end
            for i0=Is0'
                %           fprintf('  %s (%s) %s\n',signs(direc0max(i0)+1),signs(direc0(i0)+1),ttl0{i0})
                
                direc0max(i0) = mod(Direc0max(i0)+direc0(i0),2);  % used only for output
                
                fprintf('  %s (%s) %s\n',signs(direc0(i0)+1),signs(direc0max(i0)+1),ttl0{i0})
            end
        end
    end
    Signs  = signs(direc0max+1);
    Signs0 = signs(direc0+1);
    if flag.Direc0, Direc0max = Direc0; end
    
    %----------------------------------------------------------------------------------------------
    % Subtract one measure
    %----------------------------------------------------------------------------------------------
    % fprintf('\nSubtract one original measure\n')
    fprintf('\nCURRENT MEASURES: Subtracting one measure (parentheses for updated Cronbach''s alpha)\n')
    
    for j0=1:m0                                           % for each current subdomain
        if n_DS0(j0)>2                                    % if at least three measures
            Is0   = I_DS0{j0};                            % collect the measures in DS0
            CA    = fCAlpha(Dat0(:,Is0),Direc0max(Is0));  % default Cronbach's alpha
            CAs   = zeros(1,n_DS0(j0));
            for j0sub=1:n_DS0(j0)
                Is0sub     = setdiff(Is0,Is0(j0sub));
                CAs(j0sub) = fCAlpha(Dat0(:,Is0sub),Direc0max(Is0sub));
            end
            [ CAmax,Jmax ] = max(CAs);
            Imax = Is0(Jmax);
            %       fprintf('%s (%i) %5.2f (%5.2f, %3i) %s\n',u_DS0(j0,:),n_DS0(j0),CA,CAmax,Imax,ttl{Imax})
            fprintf('%s %5.2f (%5.2f) %s\n',u_DS0(j0,:),CA,CAmax,ttl0{Imax})
        end
    end
    
    %----------------------------------------------------------------------------------------------
    % CFA
    %----------------------------------------------------------------------------------------------
    % fprintf('\nCFA of original measures\n')
    fprintf('\nCURRENT MEASURES: Optimal direction versus CFA load factors\n\n')
    for j0=1:m0                                           % for each current subdomain
        if n_DS0(j0)>1                                    % if at least two measures
            Is0   = I_DS0{j0};                            % collect the measures in DS0
            CX    = nancov(Dat0(:,Is0));                  % get the covariance
            lam   = fCFA( CX );                           % get the loadings
            lam   = lam/(0.5/3);                          % scale to 1
            CA    = fCAlpha(Dat0(:,Is0),Direc0max(Is0));  % default Cronbach's alpha
            fprintf('\n%s (%.2f)\n',u_DS0(j0,:),CA)
            lam =  lam*sign(lam(1));
            for I0=1:n_DS0(j0)
                i0 = Is0(I0);
                if strcmp(Signs0(i0),'-'), lam(I0) = -lam(I0); end
                fprintf('  %s (%s) %5.2f %s\n',Signs0(i0),Signs(i0),lam(I0),ttl0{i0})
            end
        end
    end
    
    %----------------------------------------------------------------------------------------------
    % CFA/alpha with one added measure
    %----------------------------------------------------------------------------------------------
    % fprintf('\nCFA after adding one measure\n')
    fprintf('\nCURRENT MEASURES: Adding voter turnout to each subdomain (with Cronbach''s alpha before and after the addition)\n')
    ii = 42;                                              % added measure
    for j0=1:m0                                           % for each current subdomain
        Is0   = I_DS0{j0};                                % collect the measures in DS0
        datii = [Dat0(:,Is0) Dat(:,ii)];                  % append data
        CX    = nancov(datii);                            % get the covariance
        lam   = fCFA( CX );                               % get the loadings
        lam   = lam/(0.5/3);                              % scale to 1
        CA    = fCAlpha(Dat0(:,Is0),Direc0max(Is0));      % default Cronbach's alpha
        CAii0 = fCAlpha(datii,[Direc0max(Is0) 0]);        % new     Cronbach's alpha
        CAii1 = fCAlpha(datii,[Direc0max(Is0) 1]);        % new     Cronbach's alpha
        if CAii1<CAii0, Signii = '+'; else Signii = '-'; end
        fprintf('\n%s (%.2f, %.2f)\n',u_DS0(j0,:),CA,max(CAii0,CAii1))
        lam =  lam*sign(lam(1));
        for I0=1:n_DS0(j0)+1
            if I0<=n_DS0(j0)
                i0    = Is0(I0);
                if strcmp(Signs0(i0),'-'), lam(I0) = -lam(I0); end
                fprintf('  %s (%s) %5.2f %s\n',Signs0(i0),Signs(i0),lam(I0),ttl0{i0})
            else
                fprintf('  %s (%s) %5.2f %s\n',' ',Signii,lam(end),ttl{ii})
            end
        end
    end
    
    %----------------------------------------------------------------------------------------------
    % Move one measure
    %----------------------------------------------------------------------------------------------
    % fprintf('\nReplace one current measure\n')
    fprintf('\nCURRENT MEASURES: Replacing one measure with a measure from another subdomain (directionality of added measure in parenteses replaced measure on second line)\n')
    
    for j0=1:m0                                           % for each current subdomain
        if n_DS0(j0)>1                                    % if at least two measures
            Is0   = I_DS0{j0};                            % collect the measures in DS0
            CA    = fCAlpha(   Dat0(:,Is0),Direc0max(Is0));  % default Cronbach's alpha
            CAmaxs= []; Imaxs = CAmaxs; sgns = CAmaxs; i0subs = CAmaxs;
            idat  = 0;
            for i0=setdiff(1:n0,Is0)
                for i0sub=Is0'
                    idat = idat + 1;                      % increment index
                    Is0sub = setdiff(Is0,i0sub);          % remove measure from subdomain j0
                    [ CAmaxs(idat),imax,sgns(idat) ] = ...
                        fCAmax( CA, Dat0(:,Is0sub), Dat0(:,i0), 1, Direc0max(Is0sub) );
                    if imax==2                            % CA is maximum
                        Imaxs (idat) = n0+1;              % beyond last measure
                        i0subs(idat) = NaN;
                    else
                        Imaxs (idat) = i0;                % replaced measure
                        i0subs(idat) = i0sub;
                    end
                end
            end
            [ CAmax,j0submax ] = max(CAmaxs);                % find the maximum
            sgn  = sgns (j0submax);                          % get corresponding sign
            Imax = Imaxs(j0submax);
            if Imax>n0                                       % if CA is max, just list CA
                %           fprintf('%s (%i) %5.2f\n',u_DS0(j0,:),n_DS0(j0),CA)
                fprintf('\n%s %5.2f\n',u_DS0(j0,:),CA)
            else                                             % otherwise include CAmax and measure
                %           fprintf('%s (%i) %5.2f (%5.2f, %s%3i) %s\n',u_DS0(j0,:),n_DS0(j0),CA,CAmax,sgn,Imax,ttl{Imax})
                %           fprintf('                               %s\n',ttl0{i0subs(j0submax)})
                fprintf('\n%s %5.2f (%5.2f, %s) %s\n',u_DS0(j0,:),CA,CAmax,sgn,ttl0{Imax})
                fprintf('                       %s\n',ttl0{i0subs(j0submax)})
            end
        end
    end
    
    flag_print  = 1;
    I_DS0_match = fDS_match( flag_print );
    
    %----------------------------------------------------------------------------------------------
    % Replace one measure
    %----------------------------------------------------------------------------------------------
    % fprintf('\nReplace one current measure with original measure\n')
    fprintf('\nCURRENT MEASURES: Replacing one measure with an original (2014) measure from the same subdomain\n')
    for j0=1:m0                                              % for each current subdomain
        if n_DS0(j0)>1                                       % if at least two measures
            Is0   = I_DS0{j0};                               % collect the measures in DS0
            CA    = fCAlpha(Dat0(:,Is0),Direc0max(Is0));     % default Cronbach's alpha
            CAmaxs= zeros(n_DS0(j0),1); Imaxs = CAmaxs; sgns = CAmaxs;
            for j0sub=1:n_DS0(j0)                            % subtract current measure
                Is0sub = setdiff(Is0,Is0(j0sub));            % get reduced set
                [ CAmaxs(j0sub),Imaxs(j0sub),sgns(j0sub) ] = ...
                    fCAmax( CA, Dat0(:,Is0sub), Dat, I_DS0_match{j0}, Direc0max(Is0sub) );
            end
            [ CAmax,j0submax ] = max(CAmaxs);                % find the maximum
            sgn  = sgns (j0submax);                          % get corresponding sign
            Imax = Imaxs(j0submax);
            if Imax>n                                        % if CA is max, just list CA
                %           fprintf('\n%s (%i) %5.2f\n',u_DS0(j0,:),n_DS0(j0),CA)
                fprintf('\n%s %5.2f\n',u_DS0(j0,:),CA)
            else                                             % otherwise include CAmax and measure
                %           fprintf('%s (%i) %5.2f (%5.2f, %s%3i) %s\n',u_DS0(j0,:),n_DS0(j0),CA,CAmax,sgn,Imax,ttl{Imax})
                %           fprintf('                               %s\n',ttl0{Is0(j0submax)})
                fprintf('\n%s %5.2f (%5.2f, %s) %s\n',u_DS0(j0,:),CA,CAmax,sgn,ttl{Imax})
                fprintf('                        %s\n',ttl0{Is0(j0submax)})
            end
        end
    end
    
    %----------------------------------------------------------------------------------------------
    % Add one measure
    %----------------------------------------------------------------------------------------------
    % fprintf('\nAdd one original measure\n')
    fprintf('\nCURRENT MEASURES: Adding one measure with an original (2014) measure from the same subdomain\n')
    for j0=1:m0                                          % for each current subdomain
        Is0   = I_DS0{j0};                               % collect the measures in DS0
        CA    = fCAlpha(   Dat0(:,Is0),Direc0max(Is0));  % default Cronbach's alpha
        [ CAmax,Imax,sgn ] = fCAmax( CA, Dat0(:,Is0), Dat, I_DS0_match{j0}, Direc0max(Is0) );
        if Imax>n                                        % if CA is max, just list CA
            %       fprintf('%s (%i) %5.2f\n',u_DS0(j0,:),n_DS0(j0),CA)
            fprintf('\n%s %5.2f\n',u_DS0(j0,:),CA)
        else                                             % otherwise include CAmax and measure
            %       fprintf('%s (%i) %5.2f (%5.2f, %s%3i) %s\n',u_DS0(j0,:),n_DS0(j0),CA,CAmax,sgn,Imax,ttl{Imax})
            fprintf('\n%s %5.2f (%5.2f, %s) %s\n',u_DS0(j0,:),CA,CAmax,sgn,ttl{Imax})
        end
    end
    
    fprintf('\n\nNumber of Cronbach''s alpha analyses: %i\n',CAcount)
    
end

if flag.measures_corr
    %----------------------------------------------------------------------------------------------
    % Look at current measures
    %----------------------------------------------------------------------------------------------
    fprintf('\n********************************************************************************\n')
    
    for i0=1:n0                                       % for each current measure
        [ ~,Ihigh ] = sort(abs(cor0(:,i0)),'descend');% sort current correlations
        Js          = Ihigh(1:9)';                    % get indeces for top correlations
        for j0=1:n0                                   % go through all measures
            if j0~=i0 && ~ismember(j0,Js) && ismember(DS0(i0,:),DS0(j0,:),'rows')
                Js  = [Js j0];                        % append if same domain/subdomain
            end
        end
        alpmeasure(i0) = fCAlpha(Dat0(:,Js),Direc0(Js));  % Cronbach's alpha
        fprintf('\n%2i %s %s\n',round(alpmeasure(i0)*100),DS0(i0,:),ttl0{i0})
        for j=1:length(Js)                            % show the top cross-correlations
            J = Js(j);
            fprintf('  %3i %s %s\n',round(cor0(J,i0)*100),DS0(J,:),ttl0{J})
        end
    end
    
end

if flag.subtract_one
    %----------------------------------------------------------------------------------------------
    % Subtract one measure
    %----------------------------------------------------------------------------------------------
    for j0=1:m0                                       % for each DS0
        Is0 = I_DS0{j0};                              % collect the measures in DS0
        if length(Is0)>2                              % if at least 2
            for i0=1:length(Is0)                      % for each measure in DS0
                Isi = setdiff(Is0,Is0(i0));           % take it out
                alp0sub(j0,i0) = fCAlpha(...
                    Dat0(:,Isi), Direc0(Isi));        % get alpha
            end
        end
    end
    
    figure(520), clf, fig = gcf; fig.Name = 'subtract from current'; ax = gca;
    for j0=1:m0
        Is0 = I_DS0{j0};                               % collect the measures in DS0
        if length(Is0)>2                               % if at least 2
            corder = ax.ColorOrderIndex;               % get color order
            plot(Is0,alp0sub(j0,1:length(Is0)),'o-'), hold on
            ax.ColorOrderIndex = corder;
            plot(Is0,alp0(j0)+Is0*0,'--')
        end
    end
    hold off, xlabel('current measures'), ylabel('alpha')
    title('alpha after subtracting 1 current measure')
    
end

if flag.append_one
    %----------------------------------------------------------------------------------------------
    % Append one more measure
    %----------------------------------------------------------------------------------------------
    for j0=1:m0                                       % for each DS0
        for i=1:n                                     % for each original measure
            alp0add(j0,i) = fCAlpha(...
                [ Dat0(:,I_DS0{j0}) Dat(:,i)],...
                [Direc0( I_DS0{j0}) 1       ]);       % append to current, get alpha
        end
    end
    
    figure(510), fig = gcf; fig.Name = 'alp01';
    subplot(2,1,1), imagesc(alp0add), colorbar, xlabel('original measures'), ylabel('current subdomains')
    title('alpha after adding 1 original measure')
    subplot(2,1,2), plot(alp0add'), xlabel('original measures'), ylabel('alpha')
    
    figure(511), fig = gcf; fig.Name = 'alp01 one-by-one';
    for j0=1:m0                                       % for each original measure
        subplot(6,3,j0)
        plot(alp0add(j0,:))
        xlabel('original measures'), ylabel('alpha')
    end
    
end

if flag.crosscor_print
    %----------------------------------------------------------------------------------------------
    % Print highest cross-correlations
    %----------------------------------------------------------------------------------------------
    fprintf('\n********************************************************************************\n')
    
    for i0=1:n0                                       % for each current measure
        [ corhigh,Ihigh ] = sort(abs(cor(:,i0)),'descend');   % sort cross-correlations
        fprintf('\n%s %s\n',DS0(i0,:),ttl0{i0})       % print current measure labels
        for j=1:5                                     % show the top cross-correlations
            J = Ihigh(j);
            fprintf('  %3i %s %s\n',round(corhigh(j)*100),DS(J,:),ttl{J})
        end
    end
    
end

if flag.cor_plot
    %----------------------------------------------------------------------------------------------
    % Plot results
    %----------------------------------------------------------------------------------------------
    figure(500), fig = gcf; fig.Name = 'cross correlation';
    imagesc(abs(cor)), colorbar, title('cross correlation'), xlabel('current measures'), ylabel('original measures')
    
    figure(501), fig = gcf; fig.Name = 'current correlation';
    imagesc(abs(cor0)), colorbar, title('current correlation'), xlabel('current measures'), ylabel('current measures')
    
    figure(502), fig = gcf; fig.Name = 'original correlation';
    imagesc(abs(cor1)), colorbar, title('original correlation'), xlabel('original measures'), ylabel('original measures')
    
end


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
[ CAmax,IKmax ] = max([CAs(:) ; CA0]);           % take the maximum alpha
Imax = ceil(IKmax/2);                            % get the measure index (n+1 if CA is max)
Kmax = IKmax - 2*Imax + 2;                       % +/- sign for Kmax = 1,2
if Kmax==2, sgn = '-'; else sgn = ' '; end


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
        keyboard
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
fprintf('\nORIGINAL MEASURES (August 2014)\n')
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
global flag n n0 cor cor0 cor1 DS0 DS ttl0 ttl Dat0 Dat Fips
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
