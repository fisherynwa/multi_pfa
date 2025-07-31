###############################################################

ReadMe for the R code of the manuscript "Multiple multi-sample 
testing under arbitrary covariance dependency" by Vladimir Vutov and Thorsten Dickhaus.

Vladimir Vutov is the author of the R scripts.
In case of any (implementation) questions, comments, or suggestions.
Please reach out to - vkvutov@uni-bremen.de or vkvutov@gmail.com.


###############################################################

   We store the executable scripts that perform the real-data analysis of MALDI. 
   To reproduce the results illustrated in the manuscript:
   Figures (1) - (4) and Tables (5) - (8), just run the main script ‘case_study.R’.
   The dataset has been stored in the subfolder "data". 
   This data file does not cointain column names for the m/z variables, 
   but the object ‘mz_vector’ contains the relevant information. 
   For an in-depth description of this dataset, see - https://academic.oup.com/bioinformatics/article/34/7/1215/4604594.
   This dataset is publicly available at https://gitlab.informatik.uni-bremen.de/digipath/Deep_Learning_for_Tumor_Classification_in_IMS
   

 The folder structure is given as follows:

   
a) ./codes/ Here are stored the needed R scripts for
    carrying out the real-data analysis.
 - multi_pfa.R: This script contains the main function "multi_pfa()", and 
   two supporting functions: they execute 
   the multiple marginal models (MMM) and estimate the marginal covariance matrix among 
   the empirical estimates (‘mmm’ and ‘marg_cov’).
   Please refer to Section 2.3 in the manuscript. 
    
	
 - PFA.R: This function is taken from pfa::pfa.gwas().
   This script contains a function that runs 
   The principal factor approximation estimator (Section 2.5 in the manuscript).
   
   - helper_functions: Here are some helper functions. For example, "getScore" extracts the score function from each marg fit.

b) ./data/ Here is stored the MALDI data.   

c)./results/ The obtained results in terms of tables and figures will be stored
   in this subfolder.
  
  ###############################################################
  
 Description:
Program "multi_pfa()" performs the Multi-PFA approach described in the paper. 
The program considers the last category as a baseline. 

	Usage
	pfa(X, y, tval, reg, K, m)

	Arguments:
	X - A data matrix that contains the independent variables.

	y - The nominal categorical outcome. Either a vector or a data matrix that consists of the outcome of multinomial (random) counts.  

	tval - A sequence of rejection thresholds. The program accepts default parameters. 
		   Please note: Since the program aims to execute multi-sample comparisons.
		   'tval' accepts a list of sequences for each baseline-category comparison. 
		   In our example, we used a three-class outcome variable. 
		   To this end, in our script, we set two sequences for tvals.

	reg - One needs to pick which linear estimator for the common factors.
		  Here, there are two options "L2" or "L1". For more details, we refer to
		  Equations (13) - (14).

	K - the number of common factors (see Equation (10) in the manuscript).
        One needs to specify the number of common factors for each logit. However, the program accepts default values.
        Namely, one can run the software without an explicit choice.
  
		
		
      m - The subset of length m * p of the smallest (in absolute values) Z-statistics. 
	    Here, we recommend using either m = 0.9 or the default value.
	
      A general recommendation: the user can run the script first with the default values. 
       Afterward, based on the obtained results, the user can fine-tune 
       the parameters for K and tval. 
	
	A plausible way to do so, for K, is a possible inspection of the obtained eigenvalues,
	 i.e., to plot the eigenvalues. The software returns a list for each baseline-category logit.
	 For example, in our data analysis, "multi.pfa.maldi$fdp_class.1_class.3$Lamba".
	 The object lambda contains all related eigenvalues.
    	
     Regarding the rejection thresholds, our practice shows that a promising way 
      to select them by considering a huge grid of thresholds. 
	
	Value:
	
    Pvalue -  The dependency-unadjusted and two-sided p-values.
           	  Furthermore, the marginal estimates for each x-var (and for each logit pair).

    adjPvalue - The dependency-adjusted p-values.
   
    FDP -       The obtained FDP(t) results for the respective logit pair.
      
    Zvalue -    The Z-statistics, i.e., the standardized marginal estimates for this logit.
	
    Sigma  -    The estimated variance-covariance matrix among the Z-values.
	
  
  
  ###############################################################
   The code for the "case_study" script was run in the following environment:
   
R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22621)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

attached base packages:
[1] splines   stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] xtable_1.8-4  rlang_1.0.4   quantreg_5.94 SparseM_1.81  VGAM_1.1-7    gridExtra_2.3
[7] ggplot2_3.3.6 dplyr_1.0.9  

loaded via a namespace (and not attached):
 [1] pillar_1.8.0       compiler_4.2.1     tools_4.2.1        digest_0.6.29     
 [5] lifecycle_1.0.1    tibble_3.1.8       gtable_0.3.0       lattice_0.20-45   
 [9] pkgconfig_2.0.3    Matrix_1.4-1       DBI_1.1.3          cli_3.3.0         
[13] rstudioapi_0.14    withr_2.5.0        generics_0.1.3     vctrs_0.4.1       
[17] MatrixModels_0.5-0 grid_4.2.1         tidyselect_1.1.2   glue_1.6.2        
[21] R6_2.5.1           fansi_1.0.3        survival_3.3-1     farver_2.1.1      
[25] purrr_0.3.4        magrittr_2.0.3     scales_1.2.1       MASS_7.3-58       
[29] assertthat_0.2.1   colorspace_2.0-3   labeling_0.4.2     utf8_1.2.2        
[33] munsell_0.5.0     
