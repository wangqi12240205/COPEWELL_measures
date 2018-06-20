function cronbuchs_alpha
List = readtable('list_2016_12_13.csv');                        % read list table
Meas = readtable('measures2010.csv');                           % read measures
measures = table2array(Meas(:,2:end));                          % read measurs except the fips as array
SD       = readtable('subdomain2010.csv');                      % read subdomain data
subname  = unique(List.subdomain);                              % unique the name of subdomain
for sub  = 1:length(subname)                                    % for each subdomain
    
    subdomain = subname{sub};                                   % name of each subdomain
    x = subdomain;
    x(x==' ')  = [];                                            % remove the space of the name of subdomain 
    loc_domain = strcmpi(SD.Properties.VariableNames,x);        % find the location of subdomain in subdomain data
    clear ID mea_domain mea_name T var_mea covmatrix corre Id_signs var_SD
    SDdata = SD(:,loc_domain);                                  % get the data of that subdomain
    var_SD = var(table2array(SDdata),'omitnan');                % variance of subdomain
    
    ID    = List.id(strcmp(List.subdomain,subdomain));          % find the id of measures in subdomain
    signs = List.direction(strcmp(List.subdomain,subdomain));   % find the signs of measures
    
    num = length(ID);                                           % the amount of measures
    for nid = 1: num
        Id_signs{nid} = [signs{nid} num2str(ID(nid))];          % add the signs in the name of measure id
    end
    
    
    if num >1                                                   % if amount of measures is over 1
        for i = 1: num
            mea_domain(:,i) = measures(2:end,measures(1,:) ==ID(i,:));   % the data of measures in certain subdomain
            mea_name{i}     = ['x' num2str(ID(i))];                      % name of the measures
            var_mea(:,i)    = var(mea_domain(:,i),'omitnan');            % variance matrix of measures in subdomain
        end
        % covariance
        covmatrix = cov(mea_domain,'omitrows');
        % correlation
        corre = corrcov(covmatrix);
        % cronbuchs alpha
        alpha =(num./(num-1) )*(1-(sum(var_mea)./(num^2*var_SD)));
        
        % plot the figure
        figure(500+sub),
        fig=gcf; fig.Position = [100 300 1300 500];                           % figure position and size
        subplot(1,2,1)                                                        % subplot (1,2,1)
        fig = gcf; fig.Name = subdomain;                                      % figure name is subdomain name
    
        clims = [-0.5:0.25:1];                                                % color label 
        imagesc(corre,[-0.5 1]),colorbar('Ticks',clims)                       % image section with color bar
        title([subdomain ', Cronbuchs'' alpha =' num2str(round(alpha,3))],'Fontsize',18) % title
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
                %                 if mm ==1
                %                     if j == num
                %                         set(AX(j,mm),'YTick',[0 0.5])
                %                     else
                %                         if j ==1
                %                             set(AX(j,mm),'YTick',[0.5 1])
                %                         else
                %                             set(AX(j,mm),'YTick',[0.5])
                %                         end
                %                     end
                %
                %                 end
            end
        end
        
        title('Correlation Matrix','Fontsize',18)
        saveas(gca, ['figure_corr/' subdomain '.fig']);
    end
    
    
end




end