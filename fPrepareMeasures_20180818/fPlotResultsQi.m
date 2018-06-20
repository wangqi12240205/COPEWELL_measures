function fPlotResultsQi
%----------------------------------------------------------------------------------------------
% Plot results showing the selected measures and their signs and alpha
%----------------------------------------------------------------------------------------------
flag_use_composites = 1;                         % 1-use composite measures, 0-don't
if flag_use_composites
	load fPrepareMeasuresQiResults
	combination = [2 1 18 1 2 1 2 2 2 2 1];
else
	load fPrepareMeasuresQiResults_no_composites
	combination = [2 1 18 1 2 1 2 2 1 2 1];
end

for j=13:-1:8
    ttlSaved{6}{j+1} = ttlSaved{6}{j};
end
ttlSaved{6}{8} = 'C: % Leisure-time Physical Inactivity Prevalence';

sub_change	= [2 4 6 7 8 9 10 11 14 16 17];
load ttl0
flag.report = 1;                                 % 1-format for the report, 0-for analysis
datalist  = readtable('list_2016_12_13.csv');
direction = datalist.direction;
nSubdomains = length(CAsSaved);                  % number of subdomains


for mm = 1:nSubdomains
    textsubdomain = ttlSaved{mm};
    for k = 1:length(textsubdomain)
        texttemp = textsubdomain{k};
        textsplit = textscan(texttemp,'%s','Delimiter',':');
        if (textsplit{1}{1} == 'C')
            [~,jlength] = size(SignsSaved{mm});
            for j = 1:jlength
            if (SignsSaved{mm}(k,j) == 1)            
                SignsSaved{mm}(k,j) = direction{find(strcmp(textsplit{1}{2},ttl0))} == '+';
            else
                SignsSaved{mm}(k,j) = direction{find(strcmp(textsplit{1}{2},ttl0))} == '-';
            end
            end
        end
    end
    
end

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
        fig.Position(3) = 1000;
        fig.Position(4) = 800;
        
        
            
        
        subplot(2,1,1)                           % plot showing + and - and not selected
        imagesc(imagescData),hold on, ax = gca;
        cmap = [0.3, 0.3, 0.3; 1, 1, 1; 0.3, 0.3, 0.3];
        colormap(cmap);
        for ii = 1:nMeasures
            for jj = 1:nResults
                
                if (imagescData(ii,jj) == 1)
                    text(jj,ii,'+','HorizontalAlignment','center','Fontsize',10,'color','w','FontWeight','bold');
                end
                if (imagescData(ii,jj) == -1)
                    text(jj,ii,'-','HorizontalAlignment','center','Fontsize',12,'color',[1 0 0],'FontWeight','bold');
                end
            end
        end
        for ii = 1:nMeasures
            plot([0,nResults+1],[ii + 0.5,ii+ 0.5],'color',[0.7,0.7,0.7]),hold on
        end
        for jj = 1:nResults
            plot([jj + 0.5,jj+ 0.5],[0,nMeasures+1],'color',[0.7,0.7,0.7]),hold on
        end
       
        ax.YDir  = 'reverse';                    % make sure current measures are at top
        ax.YTick = 1:nMeasures;                  % label all measures
        ax.YTickLabel = ttlSaved{i};
        ax.YAxisLocation = 'right';              % put labels on the right side
        ax.XTick = 1:nResults;                  % label all measures
        ax.TickLength = [0 0];
        xlim([0,nResults]+0.5)
         for k=1:length(ax.XTick) 
            if mod(k,2)
                ax.XTickLabel{k} = sprintf('  %i',k);    % adjust the spaces before %i so that the labels line up correctly
            else
                ax.XTickLabel{k} = '';
            end
         end
        if ismember(i, sub_change)
            highlightnum = combination(sub_change == i);
            ax.XTickLabel{highlightnum} = ['\color{red} \bf' num2str(highlightnum)];
        end
        ax_height = nMeasures* 0.02;
        ax_length = nResults * 0.016;              % make the length of the plots the same
        
        ax.Position(4) = ax_height;
        ax.Position(3) = ax_length;
        ax.Position(2) = ax.Position(2) + 0.1;
        %         grid on
        
        
        subplot(2,1,2), plot(CAsSaved{i},'bo-','LineWidth',2), grid on, ax = gca; ax.XTickLabel = '';
        axis([-0.5 20.5 0 1])
        ax.Position(3) = ax_length;
        ax.Position(2) = ax.Position(2) + 0.18;
        xlim([0,nResults]+0.5)
        ylim([0.3 1])
        ylabel('                                Cronbach''s alpha')
        set(get(gca,'ylabel'),'rotation',0)
        ax.YAxisLocation = 'right';              % put labels on the right side
        ax.XTick = 0.5:nResults;                 % label all measures
        ax.XTickLabel = '';                      % no x-axis labels
        set(ax,'Fontsize',12)
        grid on
        if ~flag.report
            hold on, plot((-CAsSaved{i}+CAsSaved{i}(1))/CAsSaved{i}(1),'ro-'), hold off, ylim([0 1])
            legend('alpha','% decrease in alpha')
        end
        print(['combination/' fig.Name],'-djpeg','-r100')
    end
end