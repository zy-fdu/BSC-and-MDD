# BSC-and-MDD

Key words: brain sex continuum; major depressive disorders

Data and Code of paper:

Association between Brain Androgyny and Major Depressive Disorder in males and females: A Large-Scale Neuroimaging Analysis (Zhang et al. 2022).

The toolbox for computing the brain sex continuum could be found at https://github.com/zy-fdu/Brain-Sex-Continuum.

NOTE

1. Use of all data included in the paper was acknowledged. The YMU is the abbreviation for the Yang-Ming University dataset. The PKU is the abbreviation for the Peking University Sixth Hospital dataset.


2. The brain sex continuum was calculated by the machine learning classifier trained through the UK Biobank dataset (The Human Brain Is Best Described as Being on a Female/Male Continuum: Evidence from a Neuroimaging Connectivity Study (Zhang et al. 2021, doi: 10.1093/cercor/bhaa408), code: https://github.com/zy-fdu/Brain-Gender-Continuum). cov_YMU is a variable in table format, containing covariates of subjects from YMU dataset. The covariates including age, sex, meanFD, diagnosis, HAM-D, etc. Similarly, cov_metaMDD and cov_PKU contains covariates of subjects from metaMDD and PKU dataset respectively.

Neuroimaging data were preprocessed through Dr. Weikang Gong's resting-state fMRI preprocessing pipeline (https://github.com/weikanggong/Resting-state-fMRI-preprocessing). The pipeline descripted in the Method section were used to preprocess the data.

3.If the data and codes are used in your work, please cite the above reference, namely Association between Brain Androgyny and Major Depressive Disorder in males and females: A Large-Scale Neuroimaging Analysis (Zhang et al. 2022)

SUMMARY （Matlab2021b and R4.1.2)

step 1: Analysis of YMU data (YMU_analysis.m)

    step 1.1: Calculating functional connectivity (FC) and brain sex continuum (BSC) using YMU data.
    
    step 1.2: Calculating the classification accuracy on each sex by diagnosis group, and assess the significance through bootstrap with 10,000 repetitions.
    
    step 1.3: Two-sample t-test of the BSC between HCs and MDDs (in females and males respectively).
    
    step 1.4: Analysing the relationship between HAM-D and BSC.

step 2: Analysis of meta-MDD data (metaMDD_analysis.m)

    step 2.1: Calculating ROI-level and network-level FC, and BSC.
    
    setp 2.2: Calculating the classification accuracy on each meta-MDD site on each sex by diagnosis group.
  
    step 2.3: Calculating the percentage of subjects in each BSC group (i.e. female-like, androgynous, and male-like).
    
    step 2.4: Analysing the relationship between HAM-D and BSC on each meta-MDD site.
    
    step 2.5: Analyse the difference on FC between HCs and MDD patients in different BSC groups.
        
        step 2.5.1: Implementing ANOVAN to find the FC difference between HCs, female-like MDDs, androgynous MDDs, and male-like MDDs.
        
        step 2.5.2: Implementing t-test and bootstrap with 10,000 repetitions to calculate the FC difference between either two of the BSC groups.
        
        step 2.5.3: Calculating the difference of FC between HC and MDDs in each BSC groups.
        
        step 2.5.4: Plot the difference in each BSC group.
        
step 3: Meta-analysis (based on R 4.1.2) (BSC_MDD_meta.R)
    
    step 3.1: Input of data for meta-analysis.
    
    step 3.2: Implementing meta-for package for meta-analysis.
    
    step 3.3: Egger's test for heterogeneity.

step 4: Transcriptomic Analysis (enrichment_analysis_of_genes.r)
    
    step 4.1: calculate the weight for each ROI through meta-analysis. (FCweight_meta.R)
    
    step 4.2: AHBA data pre-processing. (S01_AHBA_datapreprocess.m)
    
    step 4.3: PLS analysis. (S02_Mediation_and_PLS_meta.m)
    
    step 4.4: Leave-One-Donor-Out analysis. (S03_leave_one_donor_out_and_PLS_meta.m)

step 5: Analysis of PKU data (BSCComparison.m)
    
    step 5.1: Calculating the BSC of PKU participants. (SVMprediction.m)

    step 5.2: Comparing patients at baseline and 8 weeks. (display_boxplot_parallel_BSC.m)
    
    step 5.3: HAM-D reduction in three BSC groups. (display_parallel_HAMD.m and display_deltaHAMD3groups.m)
    
    step 5.4: Examining the medication effect. (display_0w_8w_mean.m)
    
    step 5.5: Remission in different BSC groups. (PKU_remission.m)
