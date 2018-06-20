function fPlotResultsQi
%----------------------------------------------------------------------------------------------
% Plot results showing the selected measures and their signs and alpha
%----------------------------------------------------------------------------------------------
load fPrepareMeasuresQiResults                   % load results
nSubdomains = length(CAsSaved);                  % number of subdomains
 
for i=1:nSubdomains                              % for each measure
    nResults = length(CAsSaved{i});              % number of results
    if nResults>2                                % if at least two results
        pluses  = (~SignsSaved{i}).*includedsSaved{i};  % get the measures with + direction
        minuses = ( SignsSaved{i}).*includedsSaved{i};  % get the measures with - direction
        nMeasures   = size(pluses,1);                   % number of measures
        imagescData = zeros(size(pluses));              % set up data for imagesc
        imagescData(pluses==1) =  1;
        imagescData(minuses==1)= -1;
        
        figure(300+i), fig = gcf; fig.Name = sprintf('%s',u_DS0(i,:));
 
        subplot(3,1,1)                           % plot showing + and - and not selected
        imagesc(imagescData), ax = gca;
        ax.YDir  = 'reverse';                    % make sure current measures are at top
        ax.YTick = 1:nMeasures;                  % label all measures
        ax.YTickLabel = ttlSaved{i};
        ax.YAxisLocation = 'right';              % put labels on the right side
        ax.XTickLabel = '';                      % no x-axis labels
        ax_length = nMeasures*0.03;              % make the length of the plots the same
        ax.Position(3) = ax_length;
        grid on
        
        subplot(3,1,2:3), plot(CAsSaved{i},'bo-'), grid on, ax = gca; ax.XTickLabel = '';
        axis([-0.5 20.5 0 1])
        ax.Position(3) = ax_length;
        grid on
    end
end