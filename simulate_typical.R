axis,dataset,label_source,organism,modality,role,url,important_limitation
doublet,pbmc-ch,cell_hashing,human,scRNA-seq,locked_test,https://zenodo.org/record/4062232,Hashing detects mainly cross-sample doublets
doublet,cline-ch,cell_hashing,human,scRNA-seq,development,https://zenodo.org/record/4062232,Proliferative cell lines may differ from tissues
doublet,mkidney-ch,cell_hashing,mouse,scRNA-seq,locked_test,https://zenodo.org/record/4062232,Hashing singlets are not guaranteed homotypic-singlet truth
doublet,pbmc-2ctrl-dm,Demuxlet,human,scRNA-seq,development,https://zenodo.org/record/4062232,Genotype labels detect inter-individual mixtures
doublet,pbmc-2stim-dm,Demuxlet,human,scRNA-seq,locked_test,https://zenodo.org/record/4062232,Activation state can alter caller performance
doublet,HEK-HMEC-MULTI,MULTI-seq,human,scRNA-seq,locked_test,https://zenodo.org/record/4062232,Same-sample homotypic doublets remain incompletely labelled
mitosis,GSE121265,FUCCI_fluorescence_and_cell_time,human,STRT-seq,donor_and_plate_holdout,https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121265,FUCCI provides continuous cycle position rather than direct proof of completed cytokinesis
ploidy,GSE162959,FACS_2n_4n_DNA_content,mouse,scRNA-seq,timepoint_holdout,https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162959,4n can include G2-M diploid cells and stable tetraploid cells
