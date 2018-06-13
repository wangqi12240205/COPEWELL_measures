function fplot
clear all, close all

flag_data = 1;                     % 1-current, 2-current max, 3-final choice
flag_report = 1;                   % 0-check analysis, 1-format for report

if flag_data==1
    load dataplot_original
elseif flag_data==2
    load dataplot_currentmaxv2
elseif flag_data==3
    load dataplot_finalchoice
elseif flag_data==4
    load dataplot_finalchoice_composite
end
% 
% u_DS0 = {'CF Com','CF Eco','CF Edu','CF Foo','CF Gov','CF Hea','CF Hou', ...
%     'CF Nur','CF Pop','CF Tra','PM Eng','PM Nat','PR Pre','PV Dep','PV Ine', ...
%     'PV Vul','SC Soc'};

cmap  = colormap(parula);
nCmapUpper = 10;                                   % adjust upper part of cmap
nCmapLower = 2;                                    % adjust lower part of cmap
cmap  = [
    interp1(1:nCmapLower,cmap(1:nCmapLower,:),1:0.07:nCmapLower)
    cmap(nCmapLower+1:64-nCmapUpper,:) 
    interp1(1:nCmapUpper,cmap(64-nCmapUpper+1:end,:),1:0.25:nCmapUpper)
    ];

for j0=1:length(u_DS0)                                       % for each current subdomain
    clear ID mea_domain covmatrix corre Id_signs
    mea_domain = dataplot {j0};
    CAlast     = CA_plot{j0};
    [~,num]    = size(mea_domain);

    if num>1
        Letters = {'A','B','C','D','E','F','G','H','I','J','K'};
        
        covmatrix = cov(mea_domain,'omitrows');                               % covariance
        corre = corrcov(covmatrix);                                           % correlation
        
        figure(500+j0)                                                        % plot the figure
        fig=gcf; fig.Position = [100 300 1300 500];                           % figure position and size
        fig.Name = u_DS0{j0};
        
        subplot(1,2,1)                                                        % subplot (1,2,1)         
        clims = [-0.5:0.5:1];                                                 % color label
        imagesc(corre,[-0.5 1]), colorbar('Ticks',clims)                      % image section with color bar
        colormap(gca,cmap)

        set(gca,'xtick',[1:num]),set(gca,'xticklabel',Letters(1:num),'Fontsize',16);% xtick as A B C D...
        set(gca,'ytick',[1:num]),set(gca,'yticklabel',Letters(1:num),'Fontsize',16);% ytick as A B C D...
        
        title(sprintf('Pair-wise correlations'),'Fontsize',20)  
       
        subplot(1,2,2)                                                        % subplot(1,2,2)
        [~, AX, ~, ~, ~] = plotmatrix(mea_domain);                            % ax of figure
        set(AX,'XLim',[-0.2,1.2]);                                            % set xlim
        set(AX,'YLim',[-0.2,1.2]);                                            % set ylim
        
        set(AX(:,  1),'YTick',0.5,'Fontsize',16);                             % set xtick
        set(AX(end,:),'XTick',0.5,'Fontsize',16);                             % set ytick

        if flag_report
            title(sprintf('Cronbach''s alpha = %.2f\n',CAlast),'Fontsize',20)
        else
            title(sprintf('Cronbach''s alpha = %.2f (%.2f)\n',CAlast,calcAlpha(corre)),'Fontsize',20)
        end
        
        % set the text of correlation inside the plot matrix 
        for j = 1: num
            set(AX(j,  1),'YTickLabel',Letters{j})
            set(AX(end,j),'XTickLabel',Letters{j})
            for mm = 1: num
                if j ~=mm
                    if corre(j,mm)>0
                        co = 'black';
                        sgn= ' ';
                    else
                        co = 'red';
                        sgn= '-';
                    end
                    text(AX(j,mm),-0.1,1.0,sprintf('%s.%02i',sgn,round(abs(corre(j,mm))*100)),...
                        'color',co,'Fontsize',16)
                end
            end
        end
        if flag_data==1
    print(['figure_original/' fig.Name],'-djpeg','-r100')
elseif flag_data==2
     print(['figure_maxcurrent/' fig.Name],'-djpeg','-r100')
elseif flag_data==3
    print(['figure_finalchoice/' fig.Name],'-djpeg','-r100')
    elseif flag_data==4
    print(['figure_finalchoice_composite/' fig.Name],'-djpeg','-r100')
end
    end
    
end


function alpha = calcAlpha(corre)
%----------------------------------------------------------------------------------------------
% Calculate alpha
%----------------------------------------------------------------------------------------------
K   = size(corre,1);                             % number of measures
corrTotal = sum(corre(:));                       % sum of all correlations K^2 terms
corrOffDiagonal = corrTotal - sum(diag(corre));  % sum of the off-diagonal correlations K*(K-1) terms
alpha = (corrOffDiagonal/(K-1))/(corrTotal/K);
