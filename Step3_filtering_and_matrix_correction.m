%%
clear all

cd /Users/yuting/OneDrive/WHOI/c_comp_metabolite_analysis/CMP_exomtab_data_processing/DataprocessingRedo_March2024/Output
load CMP_exomtab_2024.06.09_D5.mat
mtabData = mtabData_D5_conc_filtered(:,1:58);
sInfo_D5 = sInfo_D5(1:58,:);

Dir = '/Users/yuting/OneDrive/WHOI/c_comp_metabolite_analysis/CMP_exomtab_data_processing/DataprocessingRedo_March2024/FiveSTDcurves';
D5_matrix_corr = readtable(string([Dir filesep 'matrix_correction_factors_D5_06092024.xlsx']));
k = contains(D5_matrix_corr.mtabNames,"cysteine dimer");
D5_matrix_corr.mtabNames(k) = regexprep(D5_matrix_corr.mtabNames(k), "cysteine dimer","cystine2");
clear k
k = contains(D5_matrix_corr.mtabNames,'3''deoxyguanosine pos');
D5_matrix_corr.mtabNames(k) = regexprep(D5_matrix_corr.mtabNames(k), '3''deoxyguanosine',"3′deoxyguanosine1");
clear k

species = {'Croco';'Ehux';'MM';'Pro';'Syn';'Tpsu'};

for i = 1:length(mtabNames_D5)
   k = strcmp(D5_matrix_corr.mtabNames,mtabNames_D5(i));
   D5_corr = D5_matrix_corr(k,:);
        for j = 1:length(species)
            s = contains(sInfo_D5.cName,species(j));
            mtabData_D5_conc_filtered_matrix(i,s)=mtabData(i,s)/D5_corr.(string(species(j)));
        end
end

if 1
    toDelete = {'inosine neg';'adenosine neg';'glyphosate pos';'glyphosate neg';...
        'pyridoxine pos';'putrescine neg';'cysteine 2 neg'};

    [c, ia, ib] = intersect(toDelete,mtabNames_D5);

    mtabData_D5_conc_filtered_matrix(ib,:)=[];
    mtabNames_D5(ib,:)=[];
    mtabData_D5_filtered(ib,:)=[];
    mtabData_D5(ib,:)=[];
    mtabData_D5_conc(ib,:)=[];
    mtabData_D5_conc_filtered(ib,:)=[];

    clear c ia ib toDelete
end

clear k s mtabData D5_corr i j Dir species 
mtabData_D5_conc_filtered_matrix(:,59:69) = mtabData_D5_conc_filtered(:,59:69);

save SkyMat_testing_3isotopes.2024.06.09_D5_matrix.mat
%%
clear all

cd /Users/yuting/OneDrive/WHOI/c_comp_metabolite_analysis/CMP_exomtab_data_processing/DataprocessingRedo_March2024/Output
load CMP_exomtab_2024.06.09_C13.mat
mtabData = mtabData_C13_conc_filtered(:,1:58);
sInfo_C13 = sInfo_C13(1:58,:);

Dir = '/Users/yuting/OneDrive/WHOI/c_comp_metabolite_analysis/CMP_exomtab_data_processing/DataprocessingRedo_March2024/FiveSTDcurves';
C13_matrix_corr = readtable(string([Dir filesep 'matrix_correction_factors_C13_04162024.xlsx']));
k = contains(C13_matrix_corr.mtabNames,"cysteine dimer");
clear k
k = contains(C13_matrix_corr.mtabNames,'3''deoxyguanosine pos');
C13_matrix_corr.mtabNames(k) = regexprep(C13_matrix_corr.mtabNames(k), '3''deoxyguanosine',"3′deoxyguanosine1");
clear k

species = {'Croco';'Ehux';'MM';'Pro';'Syn';'Tpsu'};

for i = 1:length(mtabNames_C13)
   k = strcmp(C13_matrix_corr.mtabNames,mtabNames_C13(i));
   C13_corr = C13_matrix_corr(k,:);
        for j = 1:length(species)
            s = contains(sInfo_C13.cName,species(j));
            mtabData_C13_conc_filtered_matrix(i,s)=mtabData(i,s)/C13_corr.(string(species(j)));
        end
end

if 1
        toDelete = {'inosine neg';'adenosine neg';'glyphosate pos';'glyphosate neg';...
        'pyridoxine pos';'putrescine neg';'cysteine 2 neg'};
    [c, ia, ib] = intersect(toDelete,mtabNames_C13);

    mtabData_C13_conc_filtered_matrix(ib,:)=[];
    mtabNames_C13(ib,:)=[];
    mtabData_C13_filtered(ib,:)=[];
    mtabData_C13(ib,:)=[];
    mtabData_C13_conc(ib,:)=[];
    mtabData_C13_conc_filtered(ib,:)=[];

    clear c ia ib toDelete
end

clear k s mtabData C13_corr i j Dir species 
mtabData_C13_conc_filtered_matrix(:,59:69) = mtabData_C13_conc_filtered(:,59:69);

save SkyMat_testing_3isotopes.2024.06.09_C13_matrix.mat

