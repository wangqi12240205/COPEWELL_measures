function ffinalchoicedata
%----------------------------------------------------------------------------------------------
% Prepare the data
%----------------------------------------------------------------------------------------------


CAcount = 0;
% load dataplot_current.mat
load('all_data')

load('fPrepareMeasuresQiResults_180620.mat')
	combination = [1 1 1 2 1 1 1 1 2];

sub_change  = 1:length(u_DS0);
for j =1:length(combination)
    sub_id =  sub_change(j);
    comb_id = combination(j);
    Inclu   = includedsSaved{sub_id};
    include = Inclu(:,comb_id);
    ttl_temp = ttlSaved{sub_id};
    ttl_include = ttl_temp(find(include));
    sign_temp   = SignsSaved{sub_id};
    signs_te    = sign_temp(:,comb_id);
    signs=  signs_te(find(include));
    CAs_te      = CAsSaved{sub_id};
    CAs         = CAs_te(comb_id);
    dat_temp = [];
    for k = 1:sum(include)
        texttemp = ttl_include{k};
        textsplit = textscan(texttemp,'%s','Delimiter',':');
        if (textsplit{1}{1} == 'N')
            dat_add = Dat0(:,strcmp(textsplit{1}{2},ttl0));
            if signs(k) ~= direc0(strcmp(textsplit{1}{2},ttl0))
                dat_add = 1-dat_add;
            end
            dat_temp = [dat_temp dat_add];
        else
            dat_add =Dat(:,strcmp(textsplit{1}{2},ttl));
            if signs(k) ~= direct0_current(strcmp(textsplit{1}{2},ttl))
                dat_add = 1-dat_add;
            end
            dat_temp = [dat_temp dat_add];
        end
        
    end
    dataplot{sub_id} = dat_temp;
    CA_plot{sub_id}  = CAs;
    
end

save('dataplot_finalchoice.mat','dataplot','CA_plot','u_DS0');



keyboard












