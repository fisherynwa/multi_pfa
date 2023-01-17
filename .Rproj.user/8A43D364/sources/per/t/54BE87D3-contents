################################################################################
#
#   Filename    :	case_study.R    												  
#                                                                                     
#
#   Author      :   V. Vutov      
#   Paper       :   Multiple multi-sample testing under arbitrary covariance
#                   dependency by V. Vutov and T. Dickhaus.
#
#   Purpose     :   This script executes the Multi-PFA methodology on the MALDI data
#                   that contains three cancer subtypes.
#                   The outcome variable is coded as follows: 1 - ADC, 2 - SqCC, 3 - Pancreatic 
#                   The independent values (i.e. the intensity values) are stored 
#                   in the object "mz_values", while the vector mz_vector corresponds
#                   to column names. 
#
#    R Version   :       R-4.2.1                                                                
#   
#
#   Required R packages :  dplyr, ggplot2, gridExtra, VGAM, quantreg, rlang
#
#
################################################################################

    list.of.packages <- c("dplyr", "ggplot2", 'gridExtra', 'VGAM', 'quantreg', 'rlang')
    
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    
    if(length(new.packages)) install.packages(new.packages)
    
    
    require(dplyr)
    require(ggplot2)
    require(gridExtra)
    require(VGAM)
    require(quantreg)
    require(rlang)
    
    ###########################################
    # Loading R packages, MALDI data & custom functions
    ###########################################
    
    source(file.path(".", "codes", "multi_pfa.R"))
    source(file.path(".", "codes", "helper_functions.R"))
    source(file.path(".", "codes", "PFA.R")) ## A chunk of code taken from pfa::pfa.gwas()

    load(file.path(".", "data", "maldi_data.rda"))
    
    
    
    ###########################################
    #### Applying the Multi-PFA method
    ###########################################
    
    multi.pfa.maldi <- multi_pfa(X = mz_values, y = Y, 
                                 tval = list(exp(-seq(1, 15, 0.2)), 
                                             exp(-seq(1, 17, 0.2))),
                                 reg = 'L1', K = list(3, 3))
    

    
    ###########################################
    #### Table "Main results" in the paper
    ###########################################
    
    table_results_panc_vs_adc <- multi.pfa.maldi$fdp_class.1_class.3$FDP[c(53, 55, 58, 60, 68, 71),]
                    
    table_results_panc_vs_sqcc <- multi.pfa.maldi$fdp_class.2_class.3$FDP[c(57, 59, 60, 62, 66, 79),]
        
    
   
    #######################################
    ## Figure 1 in the article           ##
    #  Examples of three mass spectrum##
    #######################################
    tiff(file = './results/Fig1.tiff', 
         width = 1800, height = 1800, units = "px",  res = 300 )
    
    par(mfrow = c(1, 3))
    
    plot(x = 520:2098, y = as.numeric(mz_values[68,]),
         type = 'l', col = 'blue', ylab = 'Intensity', xlab = 'm/z values')
    
    plot(x = 520:2098, y = as.numeric(mz_values[216,]), 
         type = 'l', col = 'blue', ylab = 'Intensity', xlab = 'm/z values')
    
    plot(x = 520:2098, y = as.numeric(mz_values[217,]), 
         type = 'l', col = 'blue', ylab = 'Intensity', xlab = 'm/z values')
    
    dev.off()
    
    ################################################
    #####  Figure 2 - Empirical Z Values        ###
    ####   Task: Pancreatic vs. ADC             ###
    ################################################
    
    df_z1_corr = data.frame(PF = multi.pfa.maldi$fdp_class.1_class.3$Zval)
    
    df_p1_corr <- data.frame(PF = multi.pfa.maldi$fdp_class.1_class.3$Pvalue$p.value)
    
    plot1 = ggplot(df_z1_corr, aes(x = PF)) + 
        geom_histogram(aes(y =..density..),
                       breaks = seq(-10, 10, by = 2), 
                       colour = "black", 
                       fill = "white") + ylab('Density') + xlab('Z-Values') +
        stat_function(fun = dnorm, args = list(mean = mean(df_z1_corr$PF), sd = sd(df_z1_corr$PF)), 
                      color = "red") 
    
    plot2 = ggplot(df_p1_corr, aes(x = PF)) + 
        geom_histogram(aes(y = ..density..),
                       breaks = seq(0, 1, by = 0.08), 
                       colour = "black", 
                       fill = "white") + xlab('P-values') + 
        theme(plot.title = element_text(hjust = 0.5)) +
        ylab('Frequency') 
    
    #pdf(file = './results/Fig2.pdf')
    tiff(file = './results/Fig2.tiff', 
        width = 1800, height = 1800, units = "px",  res = 300 )
    
    grid.arrange(plot1, plot2, ncol= 2 )
   

    dev.off()
    ################################################
    #####  Figure 3 - Empirical Z Values        ###
    ####   Task: Pancreatic vs. SqCC            ###
    ################################################
    df_z2_corr = data.frame(PF = multi.pfa.maldi$fdp_class.2_class.3$Zval)
    
    df_P2_corr <- data.frame(PF = multi.pfa.maldi$fdp_class.2_class.3$Pvalue$p.value)
    
    plot3 = ggplot(df_z2_corr, aes(x = PF)) + 
        geom_histogram(aes(y =..density..),
                       breaks = seq(-12, 12, by = 2), 
                       colour = "black", 
                       fill = "white") + ylab('Density') + xlab('Z-Values') +
        stat_function(fun = dnorm, args = list(mean = mean(df_z2_corr$PF), 
                                               sd = sd(df_z2_corr$PF)), color = "red") 
    
    plot4 = ggplot(df_P2_corr, aes(x = PF)) + 
        geom_histogram(aes(y = ..density..),
                       breaks = seq(0, 1, by = 0.08), 
                       colour = "black", 
                       fill = "white") + xlab('P-values') + 
        theme(plot.title = element_text(hjust = 0.5)) +
        #ggtitle('Task Squamous carcinoma  vs Pancreatic') +
        ylab('Frequency') 
    
    
    
    tiff(file = './results/Fig3.tiff', 
         width = 1800, height = 1800, units = "px",  res = 300 )
    
    grid.arrange(plot3, plot4, ncol= 2 )
    
    dev.off()
    
    
    ################################################
    #####  Figure 4 - Main Results              ####
    ####   Task: Pancreatic vs. ADC             ###
    ################################################
    
    number_ticks <- function(n) {function(limits) pretty(limits, n)}
    
    
    Rt_1 = multi.pfa.maldi$fdp_class.1_class.3$FDP$rejects
    t_1 = multi.pfa.maldi$fdp_class.1_class.3$FDP$t
    Vt_1 = multi.pfa.maldi$fdp_class.1_class.3$FDP$false.rejects
    FDPt_1 = multi.pfa.maldi$fdp_class.1_class.3$FDP$FDP
    FDP_3_1 = multi.pfa.maldi$fdp_class.1_class.3$FDP
    
    rt_1 = ggplot(FDP_3_1, aes(-log(t_1, base = 10), Rt_1)) + 
        geom_point()  + ylim(250, 1350) +  scale_y_continuous(breaks=number_ticks(10)) + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        ylab("R(t)") + xlab('-log(t)') 
    
    
    vt_1 = ggplot(FDP_3_1, aes(-log(t, base = 10), Vt_1)) + 
        geom_point()  + ylim(11, 1210) + scale_y_continuous(breaks=number_ticks(10)) + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        ylab("Estimated V(t)") + xlab('-log(t)')
    
    fdp_1 = ggplot(FDP_3_1, aes(-log(t, base = 10), FDPt_1)) + 
        geom_point()  + ylim(0.05, 0.9) + scale_y_continuous(breaks=number_ticks(10)) + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        ylab("Approximation of FDP(t)") + xlab('-log(t)')
    
    
    
    tiff(file = './results/Fig4.tiff', 
         width = 1800, height = 1800, units = "px",  res = 300 )
    
    grid.arrange(rt_1, vt_1, fdp_1, ncol= 3)
    
    dev.off()
    
    ################################################
    #####  Figure 5 - Main Results              ####
    ####   Task: Pancreatic vs. SqCC             ###
    ################################################
    
    
    Rt_2 = multi.pfa.maldi$fdp_class.2_class.3$FDP$rejects[1:71]
    t_2 = multi.pfa.maldi$fdp_class.2_class.3$FDP$t[1:71]
    Vt_2 = multi.pfa.maldi$fdp_class.2_class.3$FDP$false.rejects[1:71]
    FDPt_2 = multi.pfa.maldi$fdp_class.2_class.3$FDP$FDP[1:71]
    FDP_3_2 = multi.pfa.maldi$fdp_class.2_class.3$FDP[1:71,]
    
    rt_2 = ggplot(FDP_3_2, aes(-log(t_1, base = 10), Rt_2)) + 
        geom_point()  + ylim(200, 1340) +  scale_y_continuous(breaks=number_ticks(10)) + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        ylab("R(t)") + xlab('-log(t)') 
    
    
    vt_2 = ggplot(FDP_3_2, aes(-log(t, base = 10), Vt_2)) + 
        geom_point()  + ylim(11, 1210) +  scale_y_continuous(breaks=number_ticks(10)) + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        ylab("Estimated V(t)") + xlab('-log(t)')
    
    fdp_2 = ggplot(FDP_3_2, aes(-log(t, base = 10), FDPt_2)) + 
        geom_point()  + ylim(0.05, 0.9) +  scale_y_continuous(breaks=number_ticks(10)) + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        ylab("Approximation of FDP(t)") + xlab('-log(t)')
    
    
    
    #pdf(file = "./Graphs/Fig5.pdf")
    
    tiff(file = './results/Fig5.tiff', 
         width = 1800, height = 1800, units = "px",  res = 300 )
    
    
    grid.arrange(rt_2, vt_2, fdp_2, ncol= 3)
    
    dev.off()
    
    
    ################################################
    #####  Figure 6 - Zoom over a pre-defined mass
    #range with significant hypotheses for the pancreatic association.
    ################################################
    
    pdf(file = "./results/Fig6.pdf")
    sub_region = cbind(Y, mz_values[, 360:370])
    
    
    panc = sub_region[sub_region[,1] == 3,]
    
    adc = sub_region[sub_region[, 1] == 1, ]; sqcc = sub_region[sub_region[, 1] == 2,]
    
    lab = round(mz_vector[360:370,])
    
    #tiff(file = './results/Fig6.tiff', 
     #    width = 1800, height = 1800, units = "px",  res = 300 )
    
    
    par(mfrow = c(1, 2))
    
  matplot(t(panc[, 2:ncol(panc)]), type = c("l"), pch = 1, 
            col = 2, ylim = c(0, 0.006), xaxt="n", ylab = 'Intensity', xlab = 'm/z values') 
    
    par(new=TRUE)
    
    matplot(t(adc[, 2:ncol(adc)]), type = c("l"), pch = 1, 
            col = 3, ylim = c(0, 0.006), xaxt="n", ylab = 'Intensity', xlab = 'm/z values') 
    
    par(new=TRUE)
    
    matplot(t(sqcc[, 2:ncol(sqcc)]), type = c("l"), pch = 1, 
            col = 4, ylim = c(0, 0.006), xaxt="n", ylab = 'Intensity', xlab = 'm/z values') 
    
    axis(1, at=1:length(lab), labels = round(lab, 2))
    
    legend("topright", legend = c('Pancreatic', 'ADC', 'Sqcc'), col= c(2:4), pch=1) 
    
    
    sub_region_2 = cbind(Y, mz_values[, 330:340])
    
    panc_2 = sub_region_2[sub_region_2[,1] == 3,]
    
    adc_2 = sub_region_2[sub_region_2[,1] == 1,]
    
    sqcc_2 = sub_region_2[sub_region_2[,1 ] == 2,]
    
    lab_2 = round(mz_vector[330:340,])
    
    
    matplot(t(panc_2[, 2:ncol(panc_2)]), type = c("l"), pch = 1, col = 2, 
            ylim = c(0, 0.006), xaxt="n", ylab = 'Intensity') 
    
    par(new=TRUE)
    
    matplot(t(adc_2[, 2:ncol(adc_2)]), type = c("l"), pch = 1, col = 3,
            ylim = c(0, 0.006), xaxt="n", ylab = 'Intensity', xlab = 'm/z values') 
    
    par(new=TRUE)
    
    matplot(t(sqcc_2[, 2:ncol(sqcc_2)]), type = c("l"), pch = 1, col = 4, 
            ylim = c(0, 0.006), xaxt="n", ylab = 'Intensity', xlab = 'm/z values') 
    
    axis(1, at=1:length(lab_2), labels = round(lab_2, 2))
    
    legend("topright", legend = c('Pancreatic', 'ADC', 'Sqcc'), col= c(2:4), pch=1) 
    dev.off()
    
    
    ######################################################################
    
    marg_quant <- getMargQuants(y = Y, X = mz_values)
    
    ############################################
    ### Table with comparisons of FDR procedures.
    ###########################################
    
    ##########################################
    ## TABLE 7 Panc vs. ADC 
    #########################################
    
    p.val_cl1 <- multi.pfa.maldi$fdp_class.1_class.3$Pvalue$p.value
    
    mmm_cl1 <- fdr.control(p.val_cl1)
    
    marg_cl1 <- fdr.control(marg_quant$pval1)
    
    mmm_efron_cl1 <- fdr.control.efron(multi.pfa.maldi$fdp_class.1_class.3$Zval)
    
    marg_efron_cl1 <- fdr.control.efron(marg_quant$zval1)
    
    rejects_panc_vs_adc <- table_results_panc_vs_adc$rejects
    
    fdr_panc_vs_adc <- data.frame(Multi_PFA = rejects_panc_vs_adc[c(6, 4, 1)],
                                    BH = mmm_cl1$rejects_BH, 
                                    BY = mmm_cl1$rejects_BY,
                                    BH_marg = marg_cl1$rejects_BH,
                                    BY_marg = marg_cl1$rejects_BY,
                                    BH_Efron = mmm_efron_cl1$BH_Efron,
                                    BY_Efron = mmm_efron_cl1$BY_Efron,
                                    BH_Efron_marg = marg_efron_cl1$BH_Efron,
                                    BY_Efron_marg = marg_efron_cl1$BY_Efron)
    
    library(xtable)
    print(xtable(fdr_panc_vs_adc, type = "latex"), file = "./results/Table7a.tex")
    
    ##########################################
    ## TABLE 7 Panc vs. SqCC 
    #########################################
    
    p.val_cl2 <- multi.pfa.maldi$fdp_class.2_class.3$Pvalue$p.value
    
    mmm_cl2 <- fdr.control(p.val_cl2)
    
    marg_cl2 <- fdr.control(marg_quant$pval2)
    
    mmm_efron_cl2 <- fdr.control.efron(multi.pfa.maldi$fdp_class.2_class.3$Zval)
    
    marg_efron_cl2 <- fdr.control.efron(marg_quant$zval2)
    
    rejects_panc_vs_sqcc <- table_results_panc_vs_sqcc$rejects
    
    
    fdr_panc_vs_sqcc <- data.frame(Multi_PFA = rejects_panc_vs_sqcc[c(6, 5, 1)],
                                   BH = mmm_cl2$rejects_BH, 
                                   BY = mmm_cl2$rejects_BY,
                                   BH_marg = marg_cl2$rejects_BH,
                                   BY_marg = marg_cl2$rejects_BY,
                                   BH_Efron = mmm_efron_cl2$BH_Efron,
                                   BY_Efron = mmm_efron_cl2$BY_Efron,
                                   BH_Efron_marg = marg_efron_cl2$BH_Efron,
                                   BY_Efron_marg = marg_efron_cl2$BY_Efron)
    
    print(xtable(fdr_panc_vs_sqcc, type = "latex"), file = "./results/Table7b.tex")
    
    
    
    
    
    
    