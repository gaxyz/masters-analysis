computePower <- function( quantile, null, empirical ){
    null <- rev(null)
    sum(empirical > null[quantile] ) / length(empirical)
}
computeROC <- function( input_df, mcoef ){
    # compute Smax df
    max_stats<- input_df %>% 
        filter(m==mcoef) %>%
        group_by( replicate, covariance, s ) %>%
        summarise(Smax=max(hapflk))
    # Define global variables
    null_global <- max_stats %>% filter(s==0) 
    selection_global <- max_stats %>% filter(s!=0) 
    alpha <- seq(0,1,0.01)
    quantiles <- 1:length(alpha)
    # Allocate dataframe
    
    df <- data.frame(alpha) %>% as_tibble()
    # Define groupings
    covariances <- max_stats %>%
        dplyr::select(covariance) %>%
        distinct() %>%
        pull(covariance)
    
    
    for ( i in covariances ){
        
        null <- null_global %>%
            filter(covariance == i ) %>%
            pull(Smax)
        selection <- selection_global %>%
            filter(covariance == i ) %>%
            pull(Smax)
        
        qnull <- quantile(null, probs = alpha )
        power <- sapply(quantiles,
                        FUN = computePower,
                        empirical = selection,
                        null = qnull )
        
        df[i] <- power
        
    }
    
    
    df <- df %>%
        pivot_longer(cols = all_of(covariances),
                     names_to = "covariance",
                     values_to = "power"
        )
    
    return(df)
    
}

plotROC <- function(df, mcoef){
    
    auc <- df %>%
        group_by(covariance) %>% 
        summarise(AUC=sum(power)) %>%
        arrange(desc(AUC))
    df %>%
        ggplot( aes(x=alpha, y = power) ) +
        geom_path(aes(color = covariance), size=2, alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype="dashed", 
                    color="black", size=0.3) +
        theme_light()  +
        xlab("Type I Error") +
        ylab("Power") +
        xlim(c(0,1)) +
        ylim(c(0,1)) +
        labs(subtitle = paste0("m=",mcoef)) +
        scale_color_manual(values=colors) +
        annotation_custom(tableGrob(auc, rows = NULL),
                          ymax = 0.2,
                          ymin= 0.2,
                          xmax=0.8,
                          xmin=0.8
                          
        )
    
    
}

freqPWR<- function( d ,frequencies, mcoef, min_freq ,max_freq){
    # get hapflk in a certain frequency window
    selected <- frequencies %>% filter(generation==7600) %>%  
        filter(s==0.1,m==mcoef) %>%
        filter(freq<=max_freq&freq>min_freq&pop=="p5") %>% 
        select(replicate,s,m)  %>%
        distinct()%>%
        left_join(d)
    # total number of replicates for computing ROC
    
    total_replicates <- selected %>% 
        select(replicate) %>%
        distinct() %>% 
        pull() %>% 
        length()
    # get null distribution of hapflk (no selection)
    null <- d %>%
        filter(s==0, m==mcoef)
    # join two dfs together
    newdf <- full_join(selected,null)
    # compute ROC
    pW <- computeROC(newdf,mcoef) %>%
        mutate(replicates=total_replicates,
               minfrq=min_freq,
               maxfrq=max_freq) %>% mutate(m=mcoef)
    
    # return total replicates
    pW
}