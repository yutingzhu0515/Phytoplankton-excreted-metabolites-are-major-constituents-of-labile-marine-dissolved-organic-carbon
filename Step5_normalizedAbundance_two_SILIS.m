clear all

cd /Users/yuting/OneDrive/WHOI/c_comp_metabolite_analysis/DOC_TON/

DOC_data = readtable('CMP_exomtab_inventory_for_DOC_YZ_withData.2023.04.23.xlsx');
 
cd /Users/yuting/OneDrive/WHOI/c_comp_metabolite_analysis/CellCount/
count_data = readtable('CMP_phytoplankton_Exometabolites_cellCounts.xlsx');
count_data = sortrows(count_data,"species");

cd /Users/yuting/OneDrive/WHOI/c_comp_metabolite_analysis/CMP_exomtab_data_processing/

ChemicalIdentifiers = readtable('BCderv_Metabolite_ChemicalIdentifiers_with_formulas.csv');
formulas_names = ChemicalIdentifiers.Chemical_Name__From_Benzoyl_Chloride_Skyline_Transition_List_;
formulas = ChemicalIdentifiers.formula;
carbon_numbers = extractBetween(formulas,"C","H");


cd /Users/yuting/OneDrive/WHOI/c_comp_metabolite_analysis/CMP_exomtab_data_processing/DataprocessingRedo_2025

mtab_names_categories = readtable('metabolite_catogeries.xlsx');
cd Output
%this is the datafile with matrix correction and LOD filtering
load CMP_exomtab_OneMode.mat

DOC_data = sortrows(DOC_data,"SampleID");
species = {'Croco';'Ehux';'MM';'Pro';'Syn';'Tpsu'};
DOC = DOC_data.NPOC_uM_*6;
mtabNames(mtabNames=='N-acetylmuramic acid lose H2O') = 'N-acetyl-muramic acid lose H2O';


mtabData_nM = mtabData;
metabNames_for_DOC = erase( mtabNames , ' Na'); 
metabNames_for_DOC= erase( metabNames_for_DOC , 'lose H2O'); 
metabNames_for_DOC = erase( metabNames_for_DOC , 'plus H2O'); 
metabNames_for_DOC = erase( metabNames_for_DOC , 'dimer'); 

metabNames_for_DOC = regexprep(metabNames_for_DOC,'\d+$','');

metabNames_for_DOC = strrep(metabNames_for_DOC,"'","′");

for i = 1:length(mtabNames)
    name = char(metabNames_for_DOC(i));
    if name(end) == " "
        name = name(1:end-1);
    end
    if name(end) == " "
        name = name(1:end-1);
    end

    metabNames_for_DOC(i)= name ;

    k = strcmp(string(formulas_names),name);
    carbon_multiplier(1,i)= str2double(cell2mat(carbon_numbers(k)));
end

%carbon_multiplier(1,18) = carbon_multiplier(1,17);

samples = sInfo.cName(1:58);
sb = contains(samples,'B');
for i = 1:length(samples)
    if sb(i) == 1
        k = char(samples(i));
        samples(i)= cellstr(k (1:end-1));
    else
        k = char(samples(i));
        samples(i)= cellstr(k (1:end-2));
    end
end

clear k sb

data = mtabData_nM';data (isnan(data)) = 0;


[~, ~, subs] = unique (count_data.species,'first');
CellCounts = accumarray(subs, count_data.perML, [], @nanmean);
for i = 1:width(data)   
    %data_new(:,i) = accumarray(subs, data(:,i), [], @nanmean);
    %stdev_new(:,i) = accumarray(subs, data(:,i), [], @nanstd);

    for j = 1:length(species)
        k = contains(samples,species(j));
        kb = contains(samples(k),'filtrate_B');
        sub_data = data(k,i);
        sub_DOC = DOC (k);
        filtrate = sub_data(~kb);
        mediaBlank = sub_data(kb);
        filtrate_DOC = sub_DOC(~kb);
        mediaBlank_DOC = sub_DOC(kb);
        blkCorrDOC = mean(filtrate_DOC)-mean(mediaBlank_DOC);
        blkCorrmtab_data = filtrate-nanmean(mediaBlank);
        carbon_data = blkCorrmtab_data*carbon_multiplier(i);
        carbon_data_fraction = blkCorrmtab_data/blkCorrDOC/1000*carbon_multiplier(i);

        %DOCnormalized_data = blkCorrmtab_data/blkCorrDOC/1000*carbon_multiplier(i);
        cell_normalized_data =blkCorrmtab_data/CellCounts(j)/1000*1e9;
        data_new(j,i)= nanmean(blkCorrmtab_data);
        data_new_DOC(j,i)= nanmean(carbon_data);
        data_new_DOC_fraction(j,i)= nanmean(carbon_data_fraction);

        %data_new_DOC(j,i)= nanmean(DOCnormalized_data);
        data_new_cell(j,i)= nanmean(cell_normalized_data);

        stdev_new(j,i) = nanstd(blkCorrmtab_data)/nanmean(blkCorrmtab_data);
        g1 = cell2mat(repmat({'F'},length(filtrate),1));
        g2 = cell2mat(repmat({'B'},length(mediaBlank),1));
        [~,p,~] = ttest2(filtrate,mediaBlank);

        g = [g1;g2];
    
        y = [filtrate; mediaBlank];
        p1_ttest(j,i) = p;
        %p1(j,i)=anova1(y,g,'off'); 
    end
end
clear subs RSD i

blk_corr_metabData = zeros (6,length(mtabNames));
blk_corr_stdev = zeros (6,length(mtabNames));


for i = 1:length(mtabNames)

    for j = 1:length(species)
        if p1_ttest(j,i)<= 0.05 && data_new(j,i) > 0 
        p2(j,i) = 1;
        blk_corr_metabData (j,i) = data_new(j,i);
        DOC_blk_corr_metabData (j,i) = data_new_DOC(j,i);
        DOC_fraction_blk_corr_metabData (j,i) = data_new_DOC_fraction(j,i);

        cell_metabData (j,i) = data_new_cell(j,i);
        metab_RSD (j,i) = stdev_new(j,i);
        else 
        p2(j,i) = 0;
        end
    end
end

metabNames_new =  erase( mtabNames , 'pos'); 
metabNames_new =  erase( metabNames_new , 'neg'); 
[~, ~, subs] = unique (metabNames_new ,'first');
metabNames_new = unique(metabNames_new );

p2 =p2';
blk_corr_metabData = blk_corr_metabData';
DOC_blk_corr_metabData = DOC_blk_corr_metabData';
DOC_fraction_blk_corr_metabData = DOC_fraction_blk_corr_metabData';

cell_metabData = cell_metabData';
metab_RSD = metab_RSD';
blk_corr_metabData (blk_corr_metabData == 0) = NaN;
DOC_blk_corr_metabData (DOC_blk_corr_metabData == 0) = NaN;
DOC_fraction_blk_corr_metabData (DOC_fraction_blk_corr_metabData == 0) = NaN;

cell_metabData(cell_metabData==0)=NaN;
metab_RSD(metab_RSD==0)=NaN;

%S-(1,2-dicarboxyethyl) glutathione neg is removed MM it's likely below
%DOL, this  is a temporal solution before I completed the script for QC

p_new = p2;
blk_corr_metabData_ave = blk_corr_metabData;
DOC_blk_corr_metabData_ave = DOC_blk_corr_metabData;
DOC_fraction_blk_corr_metabData_ave = DOC_fraction_blk_corr_metabData;

cell_metabData_ave= cell_metabData;
metab_RSD_ave = metab_RSD;



p_sum = sum(p_new,2);


p_new(p_sum == 0,: )=[];
blk_corr_metabData_ave (p_sum == 0,: )=[];
DOC_blk_corr_metabData_ave (p_sum == 0,: )=[];
DOC_fraction_blk_corr_metabData_ave (p_sum == 0,: )=[];

cell_metabData_ave(p_sum == 0,: )=[];
metabNames_new (p_sum == 0) = [];
metab_RSD_ave(p_sum == 0,:) = [];

%need to confirm and see why putrescine and putrscine2 gives different data
if 1
    toDelete = {'serine 2';'ornithine';'cystine2';...
        'glucosamine-6-phosphate lose H20';'N-acetyl-D-glucosamine';'syringic acid';'N-acetyl-muramic acid lose H2O';...
        '4-aminobenzoic acid'};

    [c, ia, ib] = intersect(toDelete,metabNames_new);

    blk_corr_metabData_ave(ib,:)=[];
    DOC_blk_corr_metabData_ave(ib,:)=[];
    DOC_fraction_blk_corr_metabData_ave(ib,:)=[];

    cell_metabData_ave(ib,:)=[];
    metabNames_new(ib,:)=[];
    metab_RSD_ave(ib,:)=[];

    p_new (ib,:)=[];
    clear c ia ib toDelete
end

metabNames_new(metabNames_new=='glucosamine-6-phosphate') = 'glucosamine 6-phosphate';

metabNames_new(metabNames_new=='amMP ') = 'AmMP ';

metabNames_new = erase( metabNames_new , ' Na'); 
metabNames_new = erase( metabNames_new , ' lose H2O'); 
metabNames_new = erase( metabNames_new , ' plus H2O'); 
metabNames_new = regexprep(metabNames_new,'\d+$','');


species= {'Crocosphaera watsonii';'Emiliania huxleyi';...
        'Micromonas commoda';'Prochlorococcus marinus';'Synechococcus';...
        'Thalassiosira pseudonana'};

metabNames_new = strrep(metabNames_new,"'","′");

for i = 1:length(metabNames_new)  
    name = char(metabNames_new(i));
    if name(end) == " "
        name = name(1:end-1);
    end
    if name(end) == " "
        name = name(1:end-1);
    end
    metabNames_new(i)=name;
    k = strcmp(string(mtab_names_categories.mtabName),metabNames_new(i));
    clean_name(i,1) = string(mtab_names_categories.mtabName_full(k));
    mtab_catg(i,1) = string(mtab_names_categories.mtab_catg(k));
end

clear name

[SPECIES,NAME] = meshgrid(species,metabNames_new);
[~, Clean_Name]= meshgrid(species,clean_name);
[~, Mtab_Catg]= meshgrid(species, mtab_catg);

species = reshape(SPECIES,length(SPECIES)*width(SPECIES),1);
metabolites= reshape(NAME,length(NAME)*width(NAME),1);
full_names= reshape(Clean_Name,length(Clean_Name)*width(Clean_Name),1);
categories= reshape(Mtab_Catg,length(Mtab_Catg)*width(Mtab_Catg),1);
presence= reshape(p_new,length(p_new)*width(p_new),1);
abundance = reshape(blk_corr_metabData_ave,length(p_new)*width(p_new),1);
abundance_carbon = reshape(DOC_blk_corr_metabData_ave,length(p_new)*width(p_new),1);
fraction_carbon = reshape(DOC_fraction_blk_corr_metabData_ave,length(p_new)*width(p_new),1);

abundance_cell_norm = reshape(cell_metabData_ave,length(p_new)*width(p_new),1);
metab_RSD = reshape(metab_RSD_ave,length(p_new)*width(p_new),1);



PRESENCE = cell(size(presence));
PRESENCE(presence == 0) = {'N'};
PRESENCE(presence == 1) = {'Y'};

%nM_data = readtable("presenceTable.xlsx");

presence_Table= array2table(species);
presence_Table.metabolites = metabolites;
presence_Table.full_names=full_names;
presence_Table.categories=categories;
presence_Table.presence = PRESENCE;
presence_Table.metabolite_concentration= abundance;
presence_Table.carbon_concentration = abundance_carbon;
presence_Table.cell_specific_concentration = abundance_cell_norm;
presence_Table.carbon_fraction = fraction_carbon;

%presence_Table.RSD =  metab_RSD;

%presence_Table.abundance = nM_data.abundance;
%h = heatmap(presence_Table,'metabolites','species','ColorVariable','presence');
recycle on % Send to recycle bin instead of permanently deleting.
delete("normalized_abundance_LODfiltering_matrixCorr.xlsx");
writetable(presence_Table,"normalized_abundance_LODfiltering_matrixCorr.xlsx");


function y = findfirst(x)
if sum(~isnan(x))==0
    y = NaN;

else 
    x(isnan(x))=[];
    y = x(find(x,1,'first'));
end
end