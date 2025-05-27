library(tidyverse)
library(vegan)
library(dplyr)
library(naniar)
library(magrittr)
library(MatchIt)
library(tableone)
library(Matching)
library(logistf)
library(pairwiseAdonis)

#signed_transformation#--------------------
signed_log_transformation <- function(matrix_data) {
  transformed_matrix <- apply(matrix_data, c(1, 2), function(x) {
    if (x >1) {
      sign(x) * log10(abs(x))
    } else {
      -sign(x) * log10(abs(x))
    }
  })
  return(transformed_matrix)
}

#min-max normalization#------------------------
min_max_normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#automatically change string format into published format#------------
capitalize_if_needed <- function(x) {
  sapply(x, function(s) {
    s <- as.character(s)
    first_char <- substr(s, 1, 1)
    if (grepl("^[a-z]", first_char)) {
      substr(s, 1, 1) <- toupper(first_char)
    }
    return(s)
  })
}

#abundance and presence filtering#---------------------------------
filter_species <- function(data,abundance_threshold,presence_threshold) {
  #species as columns
  species_abundance <- colMeans(data)
  species_presence <- apply(data, 2, function(x) sum(x > 0) / length(x))
  df_filtered_final <- data[, species_abundance > abundance_threshold & species_presence > presence_threshold]
  return(df_filtered_final)
}

#caculate group abundance#------------------------
avg_group_abun <- function(taxa_table) {
  long_type <- pivot_longer(data = taxa_table,cols = everything(),
                            names_to = "species",values_to = "abundance")
  average_abundance <- long_type %>%
    group_by(species) %>%
    summarize(avg = mean(abundance, na.rm = TRUE))%>%arrange(desc(avg))
  return(average_abundance)
}

#generate comparison list for calculate Wilcox test#-------------------
generate_comparison_list <- function(group_vector) {
  #input a vector.  e.g:result$Group
  comparison_list <- list()
  unique_groups <- unique(group_vector)
  for (i in 1:(length(unique_groups) - 1)) {
    for (j in (i + 1):length(unique_groups)) {
      comparison_list <- append(comparison_list, list(c(unique_groups[i], unique_groups[j])))
    }
  }
  return(comparison_list)
}
#generate age cutoff for downstream analysis#------------------------------
generate_age_cutoff <- function(age_vector, change_points) {
  num_labels <- length(change_points) + 1
  if(length(change_points)>1) {
    labels <- c(paste0("<", change_points[1]), 
                paste0(change_points[-length(change_points)], "-", change_points[-1]), 
                paste0(change_points[length(change_points)],"+"))
  }else {labels <- c(paste0("<", change_points[1]), 
                     paste0(change_points[length(change_points)],"+"))}
  age_cutoff <- cut(age_vector, breaks = c(18, change_points, 100), 
                    labels = labels, include.lowest = TRUE, right = FALSE)
  return(list(cutoff=age_cutoff,labels=labels))
}

#pcoa analysis#-------------------------------------------------
pcoa_analysis <- function(data, meta_data, method = "bray", dimensions = 3) {
  set.seed(1234)
  print("Doing matrix caculating analysis,plsease keep waiting...........")
  distance_matrix <- as.matrix(vegdist(data, method = method, diag = TRUE, upper = TRUE))
  print("distance matrix generated")
  #PCoA
  pcoa_results <- cmdscale(distance_matrix, k = dimensions, eig = TRUE)
  points <- as.data.frame(pcoa_results$points)
  colnames(points) <- c("PC1", "PC2", "PC3")[1:dimensions]
  points$sample_id <- row.names(points)
  
  df_paint <- inner_join(points, meta_data, by = "sample_id")
  
  return(list(distance_mat = distance_matrix,df_paint = df_paint, eig = pcoa_results$eig))
}

#meta clean for adonis#-------------------------------------------
meta_clean_adonis <- function(meta_data, max_pct_miss = 25, 
                              remove_vars = c("subject_id", "study_condition", 
                                              "curator", "study_name","BMI",
                                              "PMID", "body_site", "sample_id", 
                                              "non_westernized","median_read_length"), 
                              na_replace_vars = c("DNA_extraction_kit")) 
{

  miss <- miss_var_summary(meta_data)

  miss_keep <- miss %>% dplyr::filter(pct_miss <= max_pct_miss)
  print(miss_keep)

  metadata_neat <- meta_data %>% dplyr::select(all_of(miss_keep$variable))
 
  row.names(metadata_neat) <- meta_data$sample_id
  
  metadata_neat <- metadata_neat %>% dplyr::select(-all_of(remove_vars))
  
  metadata_neat%<>%mutate(across(all_of(na_replace_vars), ~str_replace_na(., "Unknown")))
  
  return(list(clean_meta=metadata_neat,miss_keep=miss_keep))
}






#plot beta age cutoff boxplot#---------------------------------
plot_beta_compare <- function(group_data,var_name,reference_group,dis_max) {
  group_data[[var_name]] <- as.character(group_data[[var_name]])
  group_names <- unique(group_data[[var_name]])
  group_names <- group_names[group_names != reference_group]
  group_names <- as.list(group_names)
  group_ls <- list()
  for (i in 1:length(group_names)) {
    subset_group <- group_data%>%filter(!!sym(var_name)==reference_group|!!sym(var_name)==group_names[[i]])%>%pull(sample_id)
    dis_group <- dis_max[subset_group, subset_group]
    dis_group <- as.vector(as.dist(dis_group))
    dat <- data.frame(dis =dis_group,
                      group = rep(paste(reference_group, "vs", group_names[[i]])))
    group_ls[[i]] <- dat
  }
  paint_df <- bind_rows(group_ls)
  comp_ls <- generate_comparison_list(paint_df$group)
  
  p <- ggplot(paint_df, aes(x=group,y=dis)) + 
    geom_boxplot(aes(fill = group), width = 0.6, outlier.colour="white",outliers = F) +
    #scale_fill_manual(values = c('#e6a84b','#972b1d')) +
    #scale_y_continuous(expand = c(0,0),breaks = seq(0,1.5,0.25),limits = c(0,1.5))+
    labs(x = '', y = 'Bray-Curtis dissimilarity\n',fill='Age group')+
    geom_signif(comparisons = comp_ls, 
                map_signif_level = TRUE,
                textsize = 2.5,
                vjust = 0.3,
                step_increase = 0.2,
                test = "wilcox.test")+
    theme_classic() +
    theme(axis.title= element_text(size=14),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=10,colour="black",face="bold",hjust = 0.5)) 
  return(list(df=paint_df,plot=p))
}
#select age cutoffs and make paint data frame in disease data#--------------
select_age_cut_dis <- function(meta_data,change_points,change_merge) {
  #分组
  age_cutoff <- generate_age_cutoff(age_vector = meta_data$age,
                                    change_points = change_points)
  age_cutoff_merge <- generate_age_cutoff(age_vector = meta_data$age,
                                          change_points = change_merge)
  meta_data$age_cutoff <- age_cutoff$cutoff
  meta_data$age_cutoff_merge <- age_cutoff_merge$cutoff
  #meta_data$age_cutoff <- factor(meta_data$age_cutoff,levels = c('<26','26-46','46-65','65+'),
  #                                    ordered = TRUE)
  #meta_data$age_cutoff_merge <- factor(meta_dis_cross$age_cutoff_merge,levels = c('<46','46-65','65+'),
  #                                          ordered = TRUE)
  paint_df <- meta_data%>%filter(grepl("^(T2D|CRC|IGT|IBD|obesity)", disease))
  paint_df_main <- paint_df%>%filter(disease=="T2D"|disease=="CRC"|disease=="IBD"|
                                       disease=="IGT"|disease=="obesity")
  paint_df_main$cluster <- as.factor(paint_df_main$cluster)
  return(paint_df_main)
}
#auto plot beta age cutoff boxplot,using recycle#-------------
plot_beta_compare_dis_auto <- function(paint_df,reference_group,output,dis_mat) {
  dis_ls <- as.list(unique(paint_df$disease))
  for (i in 1:length(dis_ls)) {
    group_data_tmp <- paint_df%>%filter(disease==dis_ls[[i]])
    sampid_tmp <- group_data_tmp$sample_id
    dis_mat_tmp <- dis_mat[sampid_tmp,sampid_tmp]
    
    p2 <- plot_beta_compare(group_data = group_data_tmp,var_name = "age_cutoff_merge",
                            dis_max = dis_mat_tmp,reference_group = reference_group)
    p2 <- p2 + labs(title = dis_ls[[i]]) + theme(plot.title = element_text(hjust = 0.5))
    fname <- paste0(output,dis_ls[[i]],"_age_cutoff_merge.pdf")
    ggsave(plot = p2,filename = fname,width = 5,height = 3,units = "in",dpi = 400 )
    
    p3 <- plot_beta_compare(group_data = group_data_tmp,var_name = "age_cutoff",
                            dis_max = dis_mat_tmp,reference_group = reference_group)
    p3 <- p3 + labs(title = dis_ls[[i]]) + theme(plot.title = element_text(hjust = 0.5))
    fname1 <- paste0(output,dis_ls[[i]],"_age_cutoff.pdf")
    ggsave(plot = p3,filename = fname1,width = 5,height = 3,units = "in",dpi = 400 )
    
  }
}

#multiple groups maaslin analysis#----------------- 
multi_maaslin <- function(meta,taxa,cutoff,threshold=0.05,threds=20) {
  #taxa is the taxa data frame,must have the same rows as meta
  #cutoff is a column of the metadata,e.g "age_cutoff"
  #threshold is the threshold for FDR
  res_ls <- list()
  meta[[cutoff]] <- as.character(meta[[cutoff]])
  cutoff_ls <- as.list(unique(meta[[cutoff]]))
  pb <- txtProgressBar(min = 0, max = length(cutoff_ls), style = 3)
  start_time <- Sys.time()
  for (i in 1:length(cutoff_ls)) {
    tmp_df <- meta %>%
      mutate(tmp_lab = if_else(!!sym(cutoff) == cutoff_ls[[i]], cutoff_ls[[i]], "others"))
    
    ##Maaslin2
    output <- paste0("./Maaslin2","_res_",cutoff_ls[[i]])
    output <- as.character(output)
    ifelse(dir.exists(output),yes = print("The folder exists"),
           no = dir.create(path = output))
    fit_data = Maaslin2(
      input_data = taxa,
      input_metadata = tmp_df,
      output = output, random_effects = c("subject_id"),
      fixed_effects = c("tmp_lab","gender"),
      max_significance = threshold,correction = "BH",
      min_abundance = 0.01,min_prevalence = 0.1,cores = threds,
      reference = c("tmp_lab,others"))
    dat <- fit_data$results
    tab <- dat%>%dplyr::select(-c(N,N.not.zero,stderr,name))%>%
      filter(metadata== "tmp_lab" & qval<threshold)
    res_ls[[i]] <- tab
    # update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb) # close the bar
  res_all <- bind_rows(res_ls)
  end_time <- Sys.time()
  time_diff <- difftime(end_time, start_time, units = "mins")
  print(paste0("The analysis cost ",round(time_diff,digits = 2)," minutes"))
  print("Done Maaslin...............................")
  return(list(res_list = res_ls,res_df = res_all))
}

#calculate OR and fixed covariates#--------------------------
or_calulate <- function(dis_ls,painting_data,adjust=TRUE) {
  #painting_data usually is the combine of healthy metadata and disease metadata
  normal_facotors <- c("gender")
  vec <- unlist(strsplit(normal_facotors, ","))
  res_ls <- list()
  #country直接影响
  for (i in 1:length(dis_ls)) {
    formula <- "disease_group ~ oob_prediction"
    formula_adjust <-"disease_group~oob_prediction+gender" 
    formu <- ifelse(adjust==TRUE,yes = formula_adjust,no = formula)
    or_table <- painting_data%>%filter(disease_group=="healthy"|disease_group==dis_ls[[i]])
    or_table$disease_group <- ifelse(or_table$disease_group == "healthy", 0, 1)
    or_table[vec] <- lapply(or_table[vec], as.factor)
    mod <- glm(formula = formu,data = or_table,family = "binomial")
    or_res <- odds.ratio(mod)
    or_mage <- or_res[2,]
    if (or_mage$p < 0.001) {
      significance <- "***"
    } else if (or_mage$p < 0.01) {
      significance <- "**"
    } else if (or_mage$p < 0.05) {
      significance <- "*"
    } else {
      significance <- "ns"
    }
    or_mage$signif <- significance
    rownames(or_mage) <- dis_ls[[i]]
    res_ls[[i]] <- or_mage
  }
  or_df <- bind_rows(res_ls)
  return(or_df)
}
#firth regression#------------------------
or_calulate_firth <- function(dis_ls,painting_data,var,target,adjust=TRUE) {
  #painting_data usually is the combine of healthy metadata and disease metadata
  normal_facotors <- c("gender")
  vec <- unlist(strsplit(normal_facotors, ","))
  res_ls <- list()
  for (i in 1:length(dis_ls)) {
    formula <- paste0(target,"~",var)
    formula_adjust <- paste0(target,"~",var,"+gender")
    formu <- ifelse(adjust==TRUE,yes = formula_adjust,no = formula)
    or_table <- painting_data%>%filter(disease_group=="healthy"|disease_group==dis_ls[[i]])
    or_table$disease_group <- ifelse(or_table$disease_group == "healthy", 0, 1)
    or_table[vec] <- lapply(or_table[vec], as.factor)
    mod <- logistf(formula =formu, data = or_table)
    tmp <- summary(mod)
    coef <- tmp$coefficients[[2]]
    p_value <- tmp$prob[[2]]
    ci_up <- tmp$ci.upper[[2]]
    ci_low <- tmp$ci.lower[[2]]
    OR <- exp(coef)
    CI_up <- exp(ci_up)
    CI_low <- exp(ci_low)
    or_mage <- data.frame(OR = OR, CI_low = CI_low, 
                          CI_high = CI_up,p=p_value)
    if (or_mage$p < 0.001) {
      significance <- "***"
    } else if (or_mage$p < 0.01) {
      significance <- "**"
    } else if (or_mage$p < 0.05) {
      significance <- "*"
    } else {
      significance <- "ns"
    }
    or_mage$signif <- significance
    rownames(or_mage) <- dis_ls[[i]]
    res_ls[[i]] <- or_mage
  }
  or_df <- bind_rows(res_ls)
  return(or_df)
}


#caculate or values with random effects (glmer)#----------------------------
or_caculate_long <- function(dis_ls,painting_data,adjust=TRUE) {
  normal_facotors <- c("gender","subject_id")
  vec <- unlist(strsplit(normal_facotors, ","))
  res_ls <- list()
  for (i in 1:length(dis_ls)) {
    formula <- "disease_group ~ oob_prediction+(1|subject_id)"
    formula_adjust <-"disease_group~oob_prediction+gender+(1|subject_id)" 
    formu <- ifelse(adjust==TRUE,yes = formula_adjust,no = formula)
    or_table <- painting_data%>%filter(disease_group=="healthy"|disease_group==dis_ls[[i]])
    or_table$disease_group <- ifelse(or_table$disease_group == "healthy", 0, 1)
    or_table$disease_group <- as.factor(or_table$disease_group)
    #scale prediction value
    #or_table$oob_prediction <- scale(or_table$oob_prediction,center = TRUE,scale = TRUE)
    or_table[vec] <- lapply(or_table[vec], as.factor)
    mod <- lme4::glmer(data = or_table,formu,family = binomial(link = "logit"))
    coef_table <- summary(mod)$coefficients
    OR <- exp(coef_table[, "Estimate"][2])
    conf_int <- exp(confint(mod,method = "Wald"))
    p_value <- coef_table[2,4]
    or_mage <- data.frame(OR = OR, CI_low = conf_int[3, 1], 
                          CI_high = conf_int[3, 2],p=p_value)
    if (or_mage$p < 0.001) {
      significance <- "***"
    } else if (or_mage$p < 0.01) {
      significance <- "**"
    } else if (or_mage$p < 0.05) {
      significance <- "*"
    } else {
      significance <- "ns"
    }
    or_mage$signif <- significance
    rownames(or_mage) <- dis_ls[[i]]
    res_ls[[i]] <- or_mage
  }
  or_df <- bind_rows(res_ls)
  return(or_df)
}

#detect change points#--------------------------------------------
library(trend)
library(cpm)
analyze_cpm <- function(data, x_column, y_column, 
                        cpmType = "Kolmogorov-Smirnov",arl0=500,startup=20) {

  cpm_df_m <- data %>%
    arrange(!!sym(x_column)) %>%
    group_by(!!sym(x_column)) %>%
    summarise(mean_shannon = median(!!sym(y_column)))%>%as.data.frame()
  row.names(cpm_df_m) <- cpm_df_m[[x_column]]
  

  cpm.res <- processStream(cpm_df_m$mean_shannon, cpmType = cpmType,ARL0 = arl0,startup = startup)

  change_points <- cpm.res$changePoints +(min(cpm_df_m[[x_column]])-1)
  
  plot <- ggplot(cpm_df_m, aes(x = !!sym(x_column), y = mean_shannon)) +
    geom_line(color = "steelblue", size = 1) + xlab(x_column) + 
    ylab(paste0("Median ",y_column)) +
    geom_vline(xintercept = change_points, color = "red", size = 1)+
    theme_bw()

  print(change_points)
  
  print(length(change_points))
  
  trends <- lapply(1:(length(change_points) + 1), function(i) {
    if (i == 1) {
      trend <- data %>% filter(!!sym(x_column) < change_points[i])
      mutation_interval <- paste0("<", change_points[i])
    } else if (i == length(change_points) + 1) {
      trend <- data %>% filter(!!sym(x_column) >= change_points[i - 1])
      mutation_interval <- paste0("≥", change_points[i - 1])
    } else {
      trend <- data %>% filter(!!sym(x_column) >= change_points[i - 1] & !!sym(x_column) < change_points[i])
      mutation_interval <- paste0(change_points[i - 1], "-", change_points[i])
    }
    trend_name <- paste("trend", i, sep = "")
    mk_result <- mk.test(trend[[y_column]], continuity = TRUE)

    if (mk_result$p.value < 0.001) {
      significance <- "***"
    } else if (mk_result$p.value < 0.01) {
      significance <- "**"
    } else if (mk_result$p.value < 0.05) {
      significance <- "*"
    } else {
      significance <- "ns"
    }
    return(data.frame(Method = mk_result$method, Mutation_interval = mutation_interval,
                      Statistic = mk_result$statistic, 
                      N = mk_result$parameter, p.value = mk_result$p.value, 
                      Signif = significance,row.names = trend_name))
  })
  
  mk_results_df <- do.call(rbind, trends)
  
  return(list(change_points = change_points, mk_results = mk_results_df, 
              plot = plot))
}

#plot cpm_mage point plot
plot_cpm_mage <- function(data,x,y,change_points) {
  #data is a dataframe with age and microbiota age
  #change_points is the results of cpm analysis
  plot <- ggplot(data = data,
                 mapping = aes(x=!!sym(x),y=!!sym(y))) +
    geom_point(size=1,alpha = 0.15) + 
    labs(x="Chronological Age",y="Microbiota Age") + 
    theme_classic() +  
    theme(legend.position = c(0.82,0.2),legend.background = element_blank(),
          legend.key = element_blank(),legend.title = element_blank(),
          axis.text.x = element_text(size=8,colour="black"),
          axis.text.y = element_text(size=8,colour="black"),
          axis.title.x = element_text(size = 12,colour = "black"),
          axis.title.y = element_text(size = 12,colour = "black"))+
    geom_smooth(method = "loess",linetype=1,fullrange = FALSE,se = TRUE) +
    geom_vline(xintercept = mage_points$change_points, linetype = "dashed",
               color = "darkred", size = 1,alpha=0.5) + 
    scale_x_continuous(breaks = c(18,change_points,50,75,95))
  return(plot)
}

#pairwise.Adonis_parallel
library(pairwiseAdonis)
pairwise.adonis_prallel <- function (x, factors, sim.function = "vegdist", sim.method = "bray", 
                                     p.adjust.m = "bonferroni", reduce = NULL, perm = 999,threds=1) 
{
  co <- combn(unique(as.character(factors)), 2)
  pairs <- c()
  Df <- c()
  SumsOfSqs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  for (elem in 1:ncol(co)) {
    if (inherits(x, "dist")) {
      x1 = as.matrix(x)[factors %in% c(as.character(co[1, 
                                                       elem]), as.character(co[2, elem])), factors %in% 
                          c(as.character(co[1, elem]), as.character(co[2, 
                                                                       elem]))]
    }
    else (if (sim.function == "daisy") {
      x1 = daisy(x[factors %in% c(co[1, elem], co[2, elem]), 
      ], metric = sim.method)
    }
    else {
      x1 = vegdist(x[factors %in% c(co[1, elem], co[2, 
                                                    elem]), ], method = sim.method)
    })
    x2 = data.frame(Fac = factors[factors %in% c(co[1, elem], 
                                                 co[2, elem])])
    print("Doing adonis........................")
    ad <- adonis2(x1 ~ Fac, data = x2, permutations = perm,parallel = threds)
    print("Done adonis.........................")
    pairs <- c(pairs, paste(co[1, elem], "vs", co[2, elem]))
    Df <- c(Df, ad$Df[1])
    SumsOfSqs <- c(SumsOfSqs, ad$SumOfSqs[1])
    F.Model <- c(F.Model, ad$F[1])
    R2 <- c(R2, ad$R2[1])
    p.value <- c(p.value, ad$`Pr(>F)`[1])
  }
  p.adjusted <- p.adjust(p.value, method = p.adjust.m)
  sig = c(rep("", length(p.adjusted)))
  sig[p.adjusted <= 0.05] <- "."
  sig[p.adjusted <= 0.01] <- "*"
  sig[p.adjusted <= 0.001] <- "**"
  sig[p.adjusted <= 1e-04] <- "***"
  pairw.res <- data.frame(pairs, Df, SumsOfSqs, F.Model, R2, 
                          p.value, p.adjusted, sig)
  if (!is.null(reduce)) {
    pairw.res <- subset(pairw.res, grepl(reduce, pairs))
    pairw.res$p.adjusted <- p.adjust(pairw.res$p.value, method = p.adjust.m)
    sig = c(rep("", length(pairw.res$p.adjusted)))
    sig[pairw.res$p.adjusted <= 0.1] <- "."
    sig[pairw.res$p.adjusted <= 0.05] <- "*"
    sig[pairw.res$p.adjusted <= 0.01] <- "**"
    sig[pairw.res$p.adjusted <= 0.001] <- "***"
    pairw.res <- data.frame(pairw.res[, 1:7], sig)
  }
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
}

#build NCM model
library(Hmisc)
library(minpack.lm)
library(stats4)
NCM_build <- function(data,group_name) {
  #columns as variables, rows as samples 
  spp <- as.matrix(data)
  #Non-linear least squares，NLS
  N <- mean(apply(spp, 1, sum))
  p.m <- apply(spp, 2, mean)
  p.m <- p.m[p.m != 0]
  p <- p.m/N
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  d = 1/N
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
  m <- coef(m.fit)
  m.ci <- confint(m.fit, 'm', level=0.95)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson",return.df=TRUE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  #cleaning results
  bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
  inter.col<-rep('black',nrow(bacnlsALL))
  inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'#出现频率低于中性群落模型预测的部分
  inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'#出现频率高于中性群落模型预测的部分
  bacnlsALL$col <- inter.col
  bacnlsALL$log_p <- log10(bacnlsALL$p)
  #paint NCM result
  p1 <- ggplot(bacnlsALL, aes(x = log_p)) +
    geom_point(aes(y = freq), color = inter.col, size = 0.5) +
    geom_line(aes(y = freq.pred), color = 'blue', size = 0.5) +
    geom_line(aes(y = Lower), color = 'blue', size = 0.5, linetype = "dashed") +
    geom_line(aes(y = Upper), color = 'blue', size = 0.5, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 1.02)) +
    labs(
      x = 'Mean Relative Abundance (log10)',
      y = 'Frequency of Occurrence',title = group_name
    ) +
    theme_bw() + 
    theme(
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold", angle = 90),
      plot.title = element_text(face = "bold",hjust = 0.5),
      panel.grid = element_blank()
    ) + annotate("text", x = -Inf, y = -Inf, 
                 label = paste0(" R2"," = ", round(Rsqr, 4)), 
                 hjust = -0.1, vjust = -25) +
    annotate("text",x = -Inf,y = -Inf,label = paste0(" Nm = ", round(m * N,4)),
             hjust = -0.1, vjust = -23)
  return(list(res_df = bacnlsALL,R2 = Rsqr,m = m,figure=p1))
}


library(psych)
segment_correlation_analysis <- function(data, age_cutoff_values,age_cutoff_col, method = "spearman", adjust = "BH") {
  res_ls <- list()
  
  for (i in 1:length(age_cutoff_values)) {
    age_group <- data %>% filter(!!sym(age_cutoff_col) == age_cutoff_values[i])
    
    cor_test <- corr.test(x = age_group$age, y = age_group$oob_prediction, method = method, adjust = adjust)
    cor_res <- data.frame(Group = age_cutoff_values[i],n = cor_test$n,R = cor_test$r, P = cor_test$p, P.adj = cor_test$p.adj)
    
    res_ls[[i]] <- cor_res
    
  }
  res_all <- bind_rows(res_ls)
  return(res_all)
}

# Auto pairwise comparision
library(cocor)
pairwise_comparison <- function(results) {
  comparison_list <- list()
  unique_groups <- unique(results$Group)
  for (i in 1:(length(unique_groups) - 1)) {
    for (j in (i + 1):length(unique_groups)) {
      comparison_list <- append(comparison_list, 
                                list(paste(unique_groups[i], unique_groups[j],sep = ",")))
    }
  }
  comp_df <- data.frame()
  for (ls in 1:length(comparison_list)) {
    tmp_str <- comparison_list[[ls]]
    sp_tab <- str_split(string = tmp_str,pattern = ",",simplify = TRUE)
    comp_df <- rbind(comp_df,sp_tab)
  }
  comp_res_ls <- list()
  for (i in 1:nrow(comp_df)) {
    group1 <- comp_df$V1[i]
    group2 <- comp_df$V2[i]
    cor_res1 <- results%>%filter(Group==group1)
    cor_res2 <- results%>%filter(Group==group2)
    res <- cocor.indep.groups(r1.jk = cor_res1$R, r2.hm = cor_res2$R, 
                              n1 = cor_res1$n, n2 = cor_res2$n, data.name = c("R1", "R2"))
    tmp <- data.frame(Group=paste0(group1," vs ",group2),Z=res@fisher1925$statistic,
                      P=res@fisher1925$p.value)
    comp_res_ls <- append(comp_res_ls,list(tmp))
  }
  comp_res_df_all <- bind_rows(comp_res_ls)
  return(comp_res_df_all)
}

#roc batch analysis
library(pROC)
roc_batch <- function(dis_ls,paint_data,age_group=TRUE) {
  res_ind <- list()
  res_mage <- list()
  pic_ls <- list()
  test_res <- list()
  for (i in 1:length(dis_ls)) {
    roc_table <- paint_data%>%filter(disease_group=="healthy"|disease_group==dis_ls[[i]])
    roc_table$disease_group <- ifelse(roc_table$disease_group == "healthy", 0, 1)
    ROC<-roc(disease_group~oob_prediction+age_index,data=roc_table,
             auc=TRUE,ci=TRUE)
    p_roc <- ggroc(ROC, legacy.axes = TRUE,aes("colour"))+
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
      theme_bw()+ggtitle('ROC')+ggsci::scale_color_lancet(labels = c("Microbiota age","Age index"))+ 
      labs(title = dis_ls[[i]]) + 
      #scale_color_manual(values = c("#FFC300","#DAF7A6"))+
      theme(plot.title = element_text(hjust = 0.5)) + 
      #scale_color_discrete(labels = c("Microbiota age","Age index")) +
      guides(color = guide_legend(title = "Group")) +
      #annotate("text",x=0.85,y=0.125,label=paste("Age-AUC = ", round(ROC$age$auc,3)))+
      annotate("text",x=0.7,y=0.25,label=paste("Microbiota age-AUC = ", round(ROC$oob_prediction$auc,3)),size=3)+
      annotate("text",x=0.7,y=0.375,label=paste("Age index-AUC = ", round(ROC$age_index$auc,3)),size=3)
    pic_ls[[i]] <- p_roc
    roc_res_ind <- data.frame(AUC=as.numeric(ROC$age_index$auc),ci_low=ROC$age_index$ci[[1]],
                              ci_high=ROC$age_index$ci[[3]])
    rownames(roc_res_ind) <- dis_ls[[i]]
    roc_res_ind$group <- "Age index"
    res_ind[[i]] <- roc_res_ind
    roc_res_mage <- data.frame(AUC=as.numeric(ROC$oob_prediction$auc),ci_low=ROC$oob_prediction$ci[[1]],
                               ci_high=ROC$oob_prediction$ci[[3]])
    rownames(roc_res_mage) <- dis_ls[[i]]
    roc_res_mage$group <- "Microbiota age"
    res_mage[[i]] <- roc_res_mage
    #roc test results
    roc_test <- roc.test(ROC$age_index,ROC$oob_prediction,method="bootstrap",boot.n=1000)
    roc_test_df <- data.frame(ind_auc=as.numeric(ROC$age_index$auc),
                              mage_auc=as.numeric(ROC$oob_prediction$auc),
                              p=roc_test$p.value)
    if (roc_test_df$p < 0.001) {
      significance <- "***"
    } else if (roc_test_df$p < 0.01) {
      significance <- "**"
    } else if (roc_test_df$p < 0.05) {
      significance <- "*"
    } else {
      significance <- "ns"
    }
    roc_test_df$signif <- significance
    rownames(roc_test_df) <- dis_ls[[i]]
    test_res[[i]] <- roc_test_df
  }
  roc_ind_all <- bind_rows(res_ind)
  roc_mage_all <- bind_rows(res_mage)
  roc_test_all <- bind_rows(test_res)
  if (age_group==TRUE) {
    roc_mage_all%<>%rownames_to_column(var = "Disease")
    roc_ind_all%<>%rownames_to_column(var = "Disease")
    roc_test_all%<>%rownames_to_column(var = "Disease")
    forest_auc <- rbind(roc_mage_all,roc_ind_all)
    age_group_lab <- unique(paint_data$age_cutoff_merge)
    roc_test_all$age_group <- age_group_lab
    forest_auc$age_group <- age_group_lab
    return(list(forest_auc=forest_auc,
                auc_test=roc_test_all,pic_ls=pic_ls))
  }
  return(list(roc_ind_all=roc_ind_all,roc_mage_all=roc_mage_all,
              roc_test_all=roc_test_all,pic_ls=pic_ls))
}

#extract potential species for pglmm
extract_pglmm_data <- function(anpan_warning_path) {
  # read all lines
  lines <- readLines(anpan_warning_path, 
                     encoding = "UTF-8")  
  # find data lines
  data_line_indices <- grep("^\\s*\\d+:\\s+", lines)
  #check lines
  if (length(data_line_indices) == 0) {
    stop("Can't find potential significant species data in the file!")
  }
  
  header_line_index <- data_line_indices[1] - 2
  
  if (header_line_index <= 0) {
    header_line_index <- grep("^bug_name\\s+prop_bug\\s+prop_global", lines)
    if (length(header_line_index) == 0) {
      stop("Can't find the header")
    }
  }
  
  # extract header info
  header_line <- lines[header_line_index]
  header <- unlist(strsplit(header_line, "\\s+"))
  header <- header[header != ""]
  # extract data info
  data_lines <- lines[data_line_indices]
  data_lines <- gsub("^\\s*\\d+:\\s*", "", data_lines)
  # double clean header
  data_lines <- data_lines[!grepl("^bug_name\\s+prop_bug\\s+prop_global", data_lines)]
  data_text <- paste(data_lines, collapse = "\n")
  # create data frame
  pglmm_name <- read.table(text = data_text, header = FALSE, stringsAsFactors = FALSE)
  colnames(pglmm_name) <- header
  # clean data type
  pglmm_name$prop_bug <- as.numeric(pglmm_name$prop_bug)
  pglmm_name$prop_global <- as.numeric(pglmm_name$prop_global)
  return(pglmm_name)
}

#extract filtered anpan files for pglmm analysis
extract_pre_pglmm_files <- function(source_dir,destination_dir,name_df) {
  if (!dir.exists(destination_dir)) {
    dir.create(destination_dir, recursive = TRUE)
  }
  # Step 2: Get bug names from the data frame
  bug_names <- name_df$bug_name
  bug_names <- trimws(bug_names)
  # Step 3: List all .tsv.gz files in the source directory
  file_list <- list.files(path = source_dir, pattern ="\\.(tsv\\.gz|txt|csv|tsv)$", full.names = TRUE)
  # Step 4: Find files matching the bug names
  files_to_copy <- c()
  for (bug in bug_names) {
    # Escape special regex characters in bug names
    bug_pattern <- gsub("([.|()\\^{}+$*?]|\\[|\\]|\\/)", "\\\\\\1", bug)
    # For exact match, use '^bug_pattern$'; for partial match, just use 'bug_pattern'
    pattern <- paste0(bug_pattern,"$")
    # Find matching files
    file_basenames <- sub("\\..*$", "", basename(file_list)) 
    matching_files <- grep(pattern, file_basenames, value = TRUE)
    # Get full paths
    matching_files_full <- file_list[file_basenames %in% matching_files]
    # Add to the list
    files_to_copy <- c(files_to_copy, matching_files_full)
  }
  # Step 5: Remove duplicates
  files_to_copy <- unique(files_to_copy)
  # Step 6: Copy the files
  for (file_path in files_to_copy) {
    dest_file_path <- file.path(destination_dir, basename(file_path))
    file.copy(from = file_path, to = dest_file_path, overwrite = TRUE)
  }
  print("Files have been copy to the new place")
}

#pca function for genefamily file
pca <- function(mat) {
  centered_mat = scale(mat, scale = FALSE)
  svd_res = svd(centered_mat)
  eigs = svd_res$d^2
  tenth_max_k = sum(eigs > (eigs[1]*.1)) + 1
  message(paste0("k = ", tenth_max_k))
  d = diag(svd_res$d)[,1:tenth_max_k]
  svd_res$u %*% d
}

#batch building tree files
library(data.table)
library(ape)
library(vegan)
batch_phylo_tree <- function(data_dir, output_dir = "output", pattern = "\\.(tsv\\.gz|txt|csv|tsv)$",
                             pca = FALSE,parallel = FALSE, num_cores = NULL,dis_method="euclidean") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # get file list
  file_list <- list.files(path = data_dir, pattern = pattern, full.names = TRUE)
  process_file <- function(file,dis_method) {
    # Read files
    x <- fread(file)
    x%<>%column_to_rownames(var = "sample_id")
    x%<>%dplyr::select(-c(age,gender))
    data_matrix <- as.matrix(x)
    # transform data into presence matrix format
    data_matrix <- data_matrix*1
    if (pca == TRUE) {
      pca_df <- pca(data_matrix)
      rownames(pca_df) <- rownames(data_matrix)
      dist_matrix <- vegdist(pca_df, method = dis_method)
    } else {
      dist_matrix <- vegdist(as.data.frame(data_matrix), method = dis_method)
    }
    # build tree
    tree <- nj(dist_matrix)
    filename <- basename(file)
    filename_no_ext <- sub("\\..*$", "", filename)
    species_name <- sub("^filtered_(.*)$", "\\1", filename_no_ext)
    output_file <- paste0(output_dir, "/", species_name, "_tree.nwk")
    write.tree(tree, file = output_file)
    return(output_file)
  }
  
  if (parallel) {
    library(parallel)
    # setting cores
    if (is.null(num_cores)) {
      num_cores <- detectCores() - 1  # keep one core
    }
    
    # parallel
    if (.Platform$OS.type == "windows") {
      # Windows use parLapply
      cl <- makeCluster(num_cores)
      clusterExport(cl, varlist = c("process_file", "output_dir"), envir = environment())
      clusterEvalQ(cl, {
        library(data.table)
        library(ape)
        library(vegan)
      })
      result <- parLapply(cl, file_list, process_file)
      stopCluster(cl)
    } else {
      # Unix/Linux/Mac use mclapply
      result <- mclapply(file_list, process_file,dis_method=dis_method, mc.cores = num_cores)
    }
  } else {
    result <- lapply(file_list, process_file)
  }
  
  return(result)
}


#use raw gene abs_pres matrix to build trees,avoid errors
batch_phylo_tree_raw <- function(data_dir, labels_dir, output_dir = "output", pattern = "\\.(tsv\\.gz|txt|csv|tsv)$",
                                 pca = FALSE,rep = FALSE, parallel = FALSE, num_cores = NULL, dis_method = "euclidean") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  file_list <- list.files(path = data_dir, pattern = pattern, full.names = TRUE)
  process_file <- function(file, dis_method) {
    filename <- basename(file)
    if (rep==TRUE) {
      filename_no_ext <- sub("\\..*$", "", filename)
      filename_no_ext <- sub("filtered_","",filename_no_ext)
    }else{
      filename_no_ext <- sub("\\..*$", "", filename)
    }
    species_name <- filename_no_ext
    
    labels_file <- file.path(labels_dir, paste0("sample_labels_", species_name, ".tsv.gz"))
    
    # check if label files exist
    if (!file.exists(labels_file)) {
      warning(paste("labels file not exists:", labels_file))
      return(NULL)
    }
    
    sample_labels <- fread(labels_file)
    
    select_cols <- sample_labels[bug_well_covered == TRUE]$sample_id
    
    abd_mat <- fread(file)[ , ..select_cols] |>
      as.matrix() |>
      t()
    
    non_zero_values <- abd_mat[abd_mat != 0]
    if (length(non_zero_values) == 0) {
      warning(paste("All data in the matrix are 0,Can not caculate:", file))
      return(NULL)
    }
    half_min = min(abd_mat[abd_mat != 0]) / 2
    abd_mat[abd_mat == 0] = half_min
    labd_mat <- log(abd_mat)
    if (pca==TRUE) {
      labd_mat <- pca(labd_mat)
      rownames(labd_mat) <- rownames(abd_mat)
      dist_matrix <- vegdist(labd_mat, method = dis_method)
      tree <- dist_matrix |> ape::nj() |> ape::ladderize()
    } else{
      dist_matrix <- vegdist(labd_mat, method = dis_method)
      tree <- dist_matrix |> ape::nj() |> ape::ladderize()
    }
    output_file <- file.path(output_dir, paste0(species_name, "_tree.nwk"))
    write.tree(tree, file = output_file)
    return(output_file)
  }
  
  if (parallel) {
    library(parallel)
    # setting number of cores to use
    if (is.null(num_cores)) {
      num_cores <- detectCores() - 1  
    }
    
    # parallel
    if (.Platform$OS.type == "windows") {
      
      cl <- makeCluster(num_cores)
      clusterExport(cl, varlist = c("process_file", "output_dir", "dis_method", "labels_dir"), envir = environment())
      clusterEvalQ(cl, {
        library(data.table)
        library(ape)
        library(vegan)
        library(tibble)
      })
      result <- parLapply(cl, file_list, process_file, dis_method)
      stopCluster(cl)
    } else {
      # Unix/Linux/Mac 
      result <- mclapply(file_list, process_file, dis_method = dis_method, mc.cores = num_cores)
    }
  } else {
    result <- lapply(file_list, process_file, dis_method)
  }
  
  return(result)
}






#calculate longitude age index slope
calculate_slopes <- function(data) {
  # check necessary columns
  required_cols <- c("subject_id", "days_from_first_collection", "age_index")
  if (!all(required_cols %in% colnames(data))) {
    stop("Column subject_id,days_from_first_collection or age_index is needed!")
  }
  
  # Ensure sampling_time as numeric type,if not,turn it into numeric type data
  if (!is.numeric(data$days_from_first_collection)) {
    
    data <- data %>%
      mutate(days_from_first_collection_num = as.numeric(as.factor(days_from_first_collection)))
    sampling_time_col <- "days_from_first_collection_num"
  } else {
    sampling_time_col <- "days_from_first_collection"
  }
  
  # Group by subject_id,calculate every subject's slope
  slopes <- data %>%
    group_by(subject_id) %>%
    do({
      df = .
      # Check if age_index and sampling_time have enough points for regression (at least 2 non-empty value)
      valid_data <- df %>%
        select(age_index, all_of(sampling_time_col)) %>%
        na.omit()
      
      if (nrow(valid_data) >= 2) {
        
        model <- lm(age_index ~ ., data = valid_data)
        # extrct slopes
        slope <- coef(model)[2]
      } else {
        # non enough points
        slope <- NA
      }
      data.frame(slope = slope)
    }) %>%
    ungroup()
  
  return(slopes)
}

#Paint tree plot for single species,fan style
library(ggtreeExtra)
library(ggtree)
library(tidytree)
library(ggpubr)
library(ggnewscale)
library(anpan)
plot_pglmm_fan <- function(model_path,meta_path,tree_path,fit_path,output_dir) {
  #check the output dir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  #read files
  meta <- fread(meta_path)
  tree_file <- read.tree(file = tree_path)
  fit <- readRDS(fit_path)
  load(model_path)
  filename <- basename(tree_path)
  bug_name <- sub("\\..*$", "", filename)
  bug_name <- str_replace(bug_name,"_tree","")
  #prepare data
  meta_use <- meta%>%filter(sample_id%in%model_input$sample_id)
  p <- plot_tree_with_post(tree_file,meta,fit = fit,
                           outcome = "age_up56",labels = model_input$sample_id,
                           return_tree_df = TRUE)
  post <- p$post_df
  post_paint <- post%>%dplyr::select(label,mean)
  post_paint%<>%dplyr::rename(sample_id = label)
  tree_file$edge.length <- log(tree_file$edge.length + 1)
  tree_df <- as_tibble(tree_file)
  tree_df <- tree_df %>% 
    left_join(meta_use[,c("sample_id","age","age_up56","age_cutoff_merge",
                          "continent_merge")],
              by=c("label"="sample_id")) 
  treefile2 <- as.treedata(tree_df)
  #prepare annotation data and color panel
  age_anno <- meta_use%>%dplyr::select(c(age_cutoff_merge,sample_id,
                                         continent_merge,age_up56))
  age_anno%<>%rename(label=sample_id,age_group=age_cutoff_merge,
                     age_up = age_up56,continent = continent_merge)
  age_col <- rev(c('#C25759','#E69191','#EDB8B0','#F5DFDB','#E8EDF1','#CCE4EF',
                   '#AECFD4','#92B5CA','#599CB4'))
  age_group_col <- c("<40"="#599CB4", "40-56"="#CCE4EF", "56+"="#E69191")
  continent_col <- c("Africa"="#424242","America"="#ba4e5c",
                     "Asia"="#ceb65e","Europe"="#328eca","Oceania"="#35a276")
  age_group_logical <- c("TRUE"="#E69191","FALSE"="#CCE4EF")
  p1_tmp <- ggtree(treefile2,layout = "fan",size = 0.15, open.angle = 10) +
    xlim(-15,NA)+scale_color_gradientn(colors =age_col,name="Age")+
    geom_tippoint(mapping = aes(color=age),size=0.8,show.legend = T)
  con_lgd <- ggpubr::get_legend(p1_tmp)
  
  p1 <- ggtree(treefile2,layout = "fan",size = 0.15, open.angle = 10) +
    xlim(-15,NA)+scale_color_gradientn(colors =age_col,name="Age")+
    geom_tippoint(mapping = aes(color=age),
                  size=0.8,
                  show.legend = F)+guides(color="none")
  
  p2_tmp <- p1 +new_scale_fill()+ geom_fruit(
    data = age_anno,geom = geom_tile,
    aes(x=0,y=label,fill=factor(age_up)),
    width=5,offset=0.10,
    show.legend=T,inherit.aes=F)+
    scale_fill_manual(values = age_group_logical,name="Age above 56")
  age_lgd <- ggpubr::get_legend(p2_tmp)
  
  p2 <- p1 +new_scale_fill()+ geom_fruit(
    data = age_anno,geom = geom_tile,
    aes(x=0,y=label,fill=factor(age_up)),
    width=5,offset=0.10,
    show.legend=T,inherit.aes=F)+
    scale_fill_manual(values = age_group_logical,name="Age above 56")+
    guides(fill="none")
  
  p3_tmp <- p2 +new_scale_fill()+ geom_fruit(
    data = age_anno,geom = geom_tile,
    aes(x=0,y=label,fill=factor(continent)),
    width=5,offset=0.10,
    show.legend=T,inherit.aes=F)+
    scale_fill_manual(values = ggplot2::alpha(continent_col,0.9),name="Continent")
  con_lgd <- ggpubr::get_legend(p3_tmp)
  
  p3 <- p2 +new_scale_fill()+ geom_fruit(
    data = age_anno,geom = geom_tile,
    aes(x=0,y=label,fill=factor(continent)),
    width=5,offset=0.10,
    show.legend=T,inherit.aes=F)+
    scale_fill_manual(values = ggplot2::alpha(continent_col,0.9),name="Continent")+
    guides(fill="none")
  p4_tmp<- p3 + new_scale_fill()+
    geom_fruit(data = post_paint,geom = geom_tile,
               aes(x=0,y = sample_id,fill=mean),
               width = 5,
               offset = 0.10,
               show.legend = T,
               inherit.aes = T)+
    scale_fill_gradient2(low ="#543005",
                         mid = "#F5F5F5",
                         high ="#003C30",
                         midpoint = median(post_paint$mean),
                         name = "Posterior mean \nphylogenetic effects")
  # scale_fill_gradientn(
  #   colors = viridis::viridis(128,option = "G",direction = 1),
  #    name = "Posterior mean \nphylogenetic effects",
  #    #trans = "log10",
  #    limits = c(min(post_paint$mean),max(post_paint$mean))
  #  ) + 
  #scale_fill_viridis_c(option = "G",direction = -1,
  #                     name="Posterior mean \nphylogenetic effects",) 
  effect_lgd <- ggpubr::get_legend(p4_tmp)
  p4 <- p3 + new_scale_fill()+
    geom_fruit(data = post_paint,geom = geom_tile,
               aes(x=0,y = sample_id,fill=mean),
               width = 5,
               offset = 0.10,
               show.legend = T,
               inherit.aes = T)+
    scale_fill_gradient2(low ="#543005",
                         mid = "#F5F5F5",
                         high ="#003C30",
                         midpoint = median(post_paint$mean),
                         name = "Posterior mean \nphylogenetic effects")+ 
    guides(fill="none") + labs(title = bug_name) + 
    theme(title = element_text(face = "bold.italic",
                               color = "black",size = 12))
  
  lgd_plot <- ggarrange(age_lgd, con_lgd,effect_lgd,
                        ncol = 2,nrow = 2,widths = c(0.1,0.1),
                        heights = c(0.1,0.1))
  
  lgd_path <- paste0("./",bug_name,"_lgd.pdf")
  fan_plot_path <- paste0("./",bug_name,"_fan.pdf")
  setwd(output_dir)
  
  ggsave(plot = lgd_plot,filename = lgd_path,device = "pdf",
         units = "in",dpi = 400,width = 4,height = 4)
  ggsave(plot = p4,filename = fan_plot_path,device = "pdf",
         units = "in",dpi = 400,width = 4,height = 4)
  return(list(tree_plot=p3,legend=lgd_plot))
}

#Paint tree plot for single species,fan style,binary outcome
plot_pglmm_fan_age_group <- function(model_path,meta_path,tree_path,
                                     fit_path,output_dir) {
  #check the output dir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  #read files
  load(model_path)
  meta <- fread(meta_path)
  tree_file <- read.tree(file = tree_path)
  fit <- readRDS(fit_path)
  filename <- basename(tree_path)
  bug_name <- sub("\\..*$", "", filename)
  bug_name <- str_replace(bug_name,"_tree","")
  #prepare data
  meta_use <- meta%>%filter(sample_id%in%model_input$sample_id)
  p <- plot_tree_with_post(tree_file,meta,fit = fit,
                           outcome = "age_up56",labels = model_input$sample_id,
                           return_tree_df = TRUE)
  post <- p$post_df
  post_paint <- post%>%dplyr::select(label,mean)
  post_paint%<>%dplyr::rename(sample_id = label)
  tree_file$edge.length <- log(tree_file$edge.length + 1)
  tree_df <- as_tibble(tree_file)
  tree_df <- tree_df %>% 
    left_join(meta_use[,c("sample_id","age","age_up56","age_cutoff_merge",
                          "continent_merge")],
              by=c("label"="sample_id")) 
  treefile2 <- as.treedata(tree_df)
  #prepare annotation data and color panel
  age_anno <- meta_use%>%dplyr::select(c(age_cutoff_merge,sample_id,
                                         continent_merge,age_up56))
  age_anno%<>%rename(label=sample_id,age_group=age_cutoff_merge,
                     age_up = age_up56,continent = continent_merge)
  age_col <- rev(c('#C25759','#E69191','#EDB8B0','#F5DFDB','#E8EDF1','#CCE4EF',
                   '#AECFD4','#92B5CA','#599CB4'))
  age_group_col <- c("<40"="#599CB4", "40-56"="#CCE4EF", "56+"="#E69191")
  continent_col <- c("Africa"="#424242","America"="#ba4e5c",
                     "Asia"="#ceb65e","Europe"="#328eca","Oceania"="#35a276")
  age_group_logical <- c("TRUE"="#E69191","FALSE"="#CCE4EF")
  p1_tmp <- ggtree(treefile2,layout = "fan",size = 0.15, open.angle = 10) +
    xlim(-15,NA)+scale_color_gradientn(colors =age_col,name="Age")+
    geom_tippoint(mapping = aes(color=age),size=0.8,show.legend = T)
  con_lgd <- ggpubr::get_legend(p1_tmp)
  
  p1 <- ggtree(treefile2,layout = "fan",size = 0.15, open.angle = 10) +
    xlim(-15,NA)+scale_color_gradientn(colors =age_col,name="Age")+
    geom_tippoint(mapping = aes(color=age),
                  size=0.8,
                  show.legend = F)+guides(color="none")
  
  p2_tmp <- p1 +new_scale_fill()+ geom_fruit(
    data = age_anno,geom = geom_tile,
    aes(x=0,y=label,fill=factor(age_up)),
    width=5,offset=0.10,
    show.legend=T,inherit.aes=F)+
    scale_fill_manual(values = age_group_logical,name="Age above 56")
  age_lgd <- ggpubr::get_legend(p2_tmp)
  
  p2 <- p1 +new_scale_fill()+ geom_fruit(
    data = age_anno,geom = geom_tile,
    aes(x=0,y=label,fill=factor(age_up)),
    width=5,offset=0.10,
    show.legend=T,inherit.aes=F)+
    scale_fill_manual(values = age_group_logical,name="Age above 56")+
    guides(fill="none")
  
  p3_tmp <- p2 +new_scale_fill()+ geom_fruit(
    data = age_anno,geom = geom_tile,
    aes(x=0,y=label,fill=factor(continent)),
    width=5,offset=0.10,
    show.legend=T,inherit.aes=F)+
    scale_fill_manual(values = ggplot2::alpha(continent_col,0.9),name="Continent")
  con_lgd <- ggpubr::get_legend(p3_tmp)
  
  p3 <- p2 +new_scale_fill()+ geom_fruit(
    data = age_anno,geom = geom_tile,
    aes(x=0,y=label,fill=factor(continent)),
    width=5,offset=0.10,
    show.legend=T,inherit.aes=F)+
    scale_fill_manual(values = ggplot2::alpha(continent_col,0.9),name="Continent")+
    guides(fill="none")
  p4_tmp<- p3 + new_scale_fill()+
    geom_fruit(data = post_paint,geom = geom_tile,
               aes(x=0,y = sample_id,fill=mean),
               width = 5,
               offset = 0.10,
               show.legend = T,
               inherit.aes = T)+
    scale_fill_gradient2(low ="#543005",
                         mid = "#F5F5F5",
                         high ="#003C30",
                         midpoint = 0,
                         name = "Posterior mean \nphylogenetic effects")
  #scale_fill_viridis_c(option = "G",direction = -1,
  #                     name="Posterior mean \nphylogenetic effects",) 
  effect_lgd <- ggpubr::get_legend(p4_tmp)
  p4 <- p3 + new_scale_fill()+
    geom_fruit(data = post_paint,geom = geom_tile,
               aes(x=0,y = sample_id,fill=mean),
               width = 5,
               offset = 0.10,
               show.legend = T,
               inherit.aes = T)+
    scale_fill_gradient2(low ="#543005",
                         mid = "#F5F5F5",
                         high ="#003C30",
                         midpoint = 0,
                         name = "Posterior mean \nphylogenetic effects")+
    guides(fill="none") + labs(title = bug_name) + 
    theme(title = element_text(face = "bold.italic",
                               color = "black",size = 12))
  
  lgd_plot <- ggarrange(age_lgd, con_lgd,effect_lgd,
                        ncol = 2,nrow = 2,widths = c(0.1,0.1),
                        heights = c(0.1,0.1))
  
  lgd_path <- paste0("./",bug_name,"_lgd.pdf")
  fan_plot_path <- paste0("./",bug_name,"_fan.pdf")
  setwd(output_dir)
  
  ggsave(plot = lgd_plot,filename = lgd_path,device = "pdf",
         units = "in",dpi = 400,width = 4,height = 4)
  ggsave(plot = p4,filename = fan_plot_path,device = "pdf",
         units = "in",dpi = 400,width = 4,height = 4)
  return(list(tree_plot=p3,legend=lgd_plot))
}

#batch plot fan plots
batch_plot_fan <- function(
    species_list,
    tree_dir,
    res_dir,
    output_dir,
    meta_path = NULL,
    group=FALSE
) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (species in species_list) {
    tree_filename <- paste0(species, "_tree.nwk")
    model_filename <- paste0(tree_filename, "_inputs.RData")
    fit_filename <- paste0(tree_filename, "_pglmm_fit.RDS")  # 注意 fit_filename 包含 tree_filename
    
    # file path
    model_path <- file.path(res_dir, model_filename)
    tree_path <- file.path(tree_dir, tree_filename)
    fit_path <- file.path(res_dir, fit_filename)
    
    # check files
    missing_files <- c()
    if (!file.exists(model_path)) {
      missing_files <- c(missing_files, model_path)
    }
    if (!file.exists(tree_path)) {
      missing_files <- c(missing_files, tree_path)
    }
    if (!file.exists(fit_path)) {
      missing_files <- c(missing_files, fit_path)
    }
    
    if (length(missing_files) > 0) {
      warning(paste("Files non exists:\n", paste(missing_files, collapse = "\n")))
      next
    }
    
    tryCatch({if(group==FALSE){
      plot_pglmm_fan(
        meta_path = meta_path,
        tree_path = tree_path,
        fit_path = fit_path,
        model_path = model_path,
        output_dir = output_dir)  
    }else{
      plot_pglmm_fan_age_group(meta_path = meta_path,
                               tree_path = tree_path,
                               fit_path = fit_path,
                               model_path = model_path,
                               output_dir = output_dir)
    }
      
      cat("Plot", species, "successfully")
    }, error = function(e) {
      warning(paste("Error when processing", species, "\nError message:", e$message))
    })
  }
}

#batch extract pglmm effect for each species
batch_extract_pglmm_effect <- function(meta_path,gene_folder,
                                       labels_folder,pglmm_out_folder,
                                       taxa_list) {
  post_ls <- list()
  for (i in 1:length(taxa_list)) {
    taxa <- taxa_list[[i]]
    gene_file <- paste0(gene_folder,taxa,".tsv")
    labels_file <- paste0(labels_folder,"sample_labels_",taxa,".tsv.gz")
    fit_path <- paste0(pglmm_out_folder,taxa,"_tree.nwk_pglmm_fit.RDS")
    model_path <- paste0(pglmm_out_folder,taxa,"_tree.nwk_inputs.RData")
    
    load(model_path)
    sample_labels <- fread(labels_file)
    select_cols <- sample_labels[bug_well_covered == TRUE]$sample_id
    abd_mat <- fread(gene_file)[ , ..select_cols] |>
      as.matrix() |>
      t()
    #build trees
    non_zero_values <- abd_mat[abd_mat != 0]
    half_min = min(abd_mat[abd_mat != 0]) / 2
    abd_mat[abd_mat == 0] = half_min
    labd_mat <- log(abd_mat)
    labd_mat <- pca(labd_mat)#k=6
    rownames(labd_mat) <- rownames(abd_mat)
    dist_matrix <- vegdist(labd_mat, method = "euclidean")
    tree <- dist_matrix |> ape::nj() |> ape::ladderize()
    
    #prepare data
    meta <- fread(meta_path)
    tree_file <- tree
    fit <- readRDS(fit_path)
    #extract pglmm effect
    p <- plot_tree_with_post(tree_file,meta,fit = fit,
                             covariates = c("gender","BMI","continent"),
                             outcome = "age_up56",labels = model_input$sample_id,
                             return_tree_df = TRUE)
    post <- p$post_df
    post_paint <- post%>%dplyr::select(label,mean)
    post_paint%<>%dplyr::rename(sample_id = label)
    post_paint <- left_join(post_paint,meta[,c("sample_id","age")],by = "sample_id")
    post_paint$taxa <- taxa
    post_ls[[i]] <- post_paint
  }
  merged_post_df <- bind_rows(post_ls)
  return(merged_post_df)
}



library(clusterProfiler)
library(RColorBrewer)
library(enrichplot)
## batch gsea analysis
batch_gsea <- function(data_dir,output_dir) {
  ## map file
  pattern <- "\\.(tsv\\.gz|txt|csv|tsv)$"
  go_name_path <- "/public/data/fujingxiang/micro_aging/mapping_v201901b/map_go_name.txt.gz"
  map_go <- "/public/data/fujingxiang/micro_aging/mapping_v201901b/map_go_uniref90_long.txt"
  go_match_uniref <- fread(map_go,header = F, sep = "\t")
  colnames(go_match_uniref) <- c("go","gene")
  #go name file
  go_name <- fread(go_name_path,header = F)
  colnames(go_name) <- c("go","term")
  
  ## data dir are those files in the model_stats dir of anpan results 
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  file_list <- list.files(path = data_dir, pattern = pattern, full.names = TRUE)
  for (file in file_list) {
    filename <- basename(file)
    filename_no_ext <- sub("\\..*$", "", filename)
    bugName <- sub("_gene_terms", "\\1", filename_no_ext)
    
    ## anpan res
    anpan_res <- fread(file)
    anpan_res <- left_join(anpan_res,go_match_uniref,by = "gene")
    anpan_res$go <- str_replace_na(anpan_res$go,replacement = "UNMAPPED")
    #filtered unmapped genes
    anpan_go_neat <- anpan_res%>%filter(go!="UNMAPPED")
    anpan_go_list <- anpan_go_neat%>%group_by(gene)%>%
      summarize(statistic=mean(statistic))
    #gene rank list
    gene_list <- anpan_go_list$statistic
    names(gene_list) <- anpan_go_list$gene
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    #select go match
    go_term_gene_df <- go_match_uniref[gene%chin%anpan_res$gene,]
    go_term_name_df <- go_name[go%chin%go_term_gene_df$go,]
    
    #GSEA
    gsea_res <- GSEA(geneList = gene_list,by = "fgsea",
                     seed = TRUE,TERM2GENE = go_term_gene_df,
                     TERM2NAME = go_term_name_df,
                     pAdjustMethod = "BH",pvalueCutoff = 0.25,
                     minGSSize = 3,maxGSSize = Inf,nPermSimple = 10000)
    setwd(output_dir)
    result <- data.frame(gsea_res)
    if (nrow(result) == 0) {
      log_message <- paste0("Species", bugName, "no term enriched under specific pvalueCutoff")
      log_file <- file.path(output_dir, paste0("gsea_", bugName, "_log.txt"))
      write(log_message, file = log_file)
      next
    }
    res_dir <- paste0("./",bugName)
    if (!dir.exists(res_dir)) {
      dir.create(res_dir, recursive = TRUE)
    }
    setwd(res_dir)
    result$species <- bugName
    df_name <- paste0("gsea_",bugName,".csv")
    write.csv(result,file = df_name,row.names = FALSE)
    
    for (i in 1:nrow(result)) {
      des_info <- result$Description[i]
      title <- paste0(bugName,"\n",des_info)
      p <- gseaplot2(gsea_res,geneSetID = rownames(result)[i],
                     title = title,
                     base_size = 7,
                     rel_heights = c(1, 0.2, 0.4),
                     subplots = 1:3,
                     pvalue_table = FALSE,
                     ES_geom = "line")
      
      plot_name <- paste0(bugName,"_gsea_",des_info,".pdf")
      ggsave(filename = plot_name,plot = p,width = 4,height = 3,
             units = "in",device = "pdf",dpi = 400)
    }
    
  }
}

#Fixed webr package PieDonut function for custom color settings
library(grid)
library(ggforce)
library(moonBook)
library(webr)
PieDonut2<-function(data,mapping,start=getOption("PieDonut.start",0),
                    addPieLabel=TRUE,addDonutLabel=TRUE,showRatioDonut=TRUE,
                    showRatioPie=TRUE,ratioByGroup=TRUE,
                    showRatioThreshold=getOption("PieDonut.showRatioThreshold",0.02),
                    labelposition=getOption("PieDonut.labelposition",2),
                    labelpositionThreshold=0.1,r0=getOption("PieDonut.r0",0.3),
                    r1=getOption("PieDonut.r1",1),r2=getOption("PieDonut.r2",1.2),
                    explode=NULL,selected=NULL,explodePos=0.1,
                    color="white",pieAlpha=0.8,donutAlpha=1,maxx=NULL,
                    showPieName=TRUE,showDonutName=FALSE,title=NULL,
                    pieLabelSize=4,donutLabelSize=3,titlesize=5,explodePie=TRUE,
                    explodeDonut=FALSE,use.label=TRUE,use.labels=TRUE,
                    family=getOption("PieDonut.family",""),man_color)
{
  (cols=colnames(data))
  if (use.labels)
    data=addLabelDf(data,mapping)
  count<-NULL
  if ("count"%in%names(mapping))
    count<-getMapping(mapping,"count")
  count
  pies<-donuts<-NULL
  (pies=getMapping(mapping,"pies"))
  if (is.null(pies))
    (pies=getMapping(mapping,"pie"))
  if(is.null(pies))
    (pies=getMapping(mapping,"x"))
  (donuts=getMapping(mapping,"donuts"))
  if (is.null(donuts))
    (donuts=getMapping(mapping,"donut"))
  if (is.null(donuts))
    (donuts=getMapping(mapping,"y"))
  if (!is.null(count)) {
    df<-data%>%group_by(.data[[pies]])%>%dplyr::summarize(Freq=sum(.data[[count]]))
    df
  }
  else {
    df=data.frame(table(data[[pies]]))
  }
  colnames(df)[1]=pies
  df$end=cumsum(df$Freq)
  df$start=dplyr::lag(df$end)
  df$start[1]=0
  total=sum(df$Freq)
  df$start1=df$start*2*pi/total
  df$end1=df$end*2*pi/total
  df$start1=df$start1+start
  df$end1=df$end1+start
  df$focus=0
  if (explodePie)
    df$focus[explode]=explodePos
  df$mid=(df$start1+df$end1)/2
  df$x=ifelse(df$focus==0,0,df$focus*sin(df$mid))
  df$y=ifelse(df$focus==0,0,df$focus*cos(df$mid))
  df$label=df[[pies]]
  df$ratio=df$Freq/sum(df$Freq)
  if (showRatioPie) {
    df$label=ifelse(df$ratio>=showRatioThreshold,
                    paste0(df$label,"\n(",scales::percent(df$ratio),")"),
                    as.character(df$label))
  }
  df$labelx=(r0+r1)/2*sin(df$mid)+df$x
  df$labely=(r0+r1)/2*cos(df$mid)+df$y
  if (!is.factor(df[[pies]]))
    df[[pies]]<-factor(df[[pies]])
  df
  mainCol=man_color
  df$radius=r1
  df$radius[df$focus!=0]=df$radius[df$focus!=0]+df$focus[df$focus!=0]
  df$hjust=ifelse((df$mid%%(2*pi))>pi,1,0)
  df$vjust=ifelse(((df$mid%%(2*pi))<(pi/2))|(df$mid%%(2*pi)>(pi*3/2)),0,1)
  df$segx=df$radius*sin(df$mid)
  df$segy=df$radius*cos(df$mid)
  df$segxend=(df$radius+0.05)*sin(df$mid)
  df$segyend=(df$radius+0.05)*cos(df$mid)
  df
  if (!is.null(donuts)) {
    subColor=makeSubColor(mainCol,no=length(unique(data[[donuts]])))
    subColor
    data
    if(!is.null(count)) {
      df3<-as.data.frame(data[c(donuts,pies,count)])
      colnames(df3)=c("donut","pie","Freq")
      df3
      df3<-eval(parse(text="complete(df3,donut,pie)"))
      df3$Freq[is.na(df3$Freq)]=0
      if (!is.factor(df3[[1]]))
        df3[[1]]=factor(df3[[1]])
      if (!is.factor(df3[[2]]))
        df3[[2]]=factor(df3[[2]])
      df3<-df3%>%arrange(.data$pie,.data$donut)
      a<-df3%>%spread(.data$pie,value=.data$Freq)
      a=as.data.frame(a)
      a
      rownames(a)=a[[1]]
      a=a[-1]
      a
      colnames(df3)[1:2]=c(donuts,pies)
    }
    else{
      df3=data.frame(table(data[[donuts]],data[[pies]]),
                     stringsAsFactors=FALSE)
      colnames(df3)[1:2]=c(donuts,pies)
      a=table(data[[donuts]],data[[pies]])
      a
    }
    a
    df3
    df3$group=rep(colSums(a),each=nrow(a))
    df3$pie=rep(1:ncol(a),each=nrow(a))
    total=sum(df3$Freq)
    total
    df3$ratio1=df3$Freq/total
    df3
    if (ratioByGroup){
      df3$ratio=scales::percent(df3$Freq/df3$group)
    }
    else{
      df3$ratio<-scales::percent(df3$ratio1)
    }
    df3$end=cumsum(df3$Freq)
    df3
    df3$start=dplyr::lag(df3$end)
    df3$start[1]=0
    df3$start1=df3$start*2*pi/total
    df3$end1=df3$end*2*pi/total
    df3$start1=df3$start1+start
    df3$end1=df3$end1+start
    df3$mid=(df3$start1+df3$end1)/2
    df3$focus=0
    if (!is.null(selected)){
      df3$focus[selected]=explodePos
    }
    else if (!is.null(explode)) {
      selected=c()
      for (i in 1:length(explode)) {
        start=1+nrow(a)*(explode[i]-1)
        selected=c(selected,start:(start+nrow(a)-1))
      }
      selected
      df3$focus[selected]=explodePos
    }
    df3
    df3$x=0
    df3$y=0
    df
    if (!is.null(explode)) {
      explode
      for (i in 1:length(explode)) {
        xpos=df$focus[explode[i]]*sin(df$mid[explode[i]])
        ypos=df$focus[explode[i]]*cos(df$mid[explode[i]])
        df3$x[df3$pie==explode[i]]=xpos
        df3$y[df3$pie==explode[i]]=ypos
      }
    }
    df3$no=1:nrow(df3)
    df3$label=df3[[donuts]]
    if (showRatioDonut) {
      if (max(nchar(levels(df3$label))) <= 2)
        df3$label=paste0(df3$label,"(",df3$ratio,
                         ")")
      else df3$label=paste0(df3$label,"\n(",df3$ratio,")")
    }
    df3$label[df3$ratio1==0]=""
    df3$label[df3$ratio1<showRatioThreshold]=""
    df3$hjust=ifelse((df3$mid%%(2*pi))>pi,1,0)
    df3$vjust=ifelse(((df3$mid%%(2*pi))<(pi/2))|
                       (df3$mid%%(2*pi)>(pi*3/2)),0,1)
    df3$no=factor(df3$no)
    df3
    labelposition
    if (labelposition>0) {
      df3$radius=r2
      if (explodeDonut)
        df3$radius[df3$focus!=0]=df3$radius[df3$focus!=0]+df3$focus[df3$focus!=0]
      df3$segx=df3$radius*sin(df3$mid)+df3$x
      df3$segy=df3$radius*cos(df3$mid)+df3$y
      df3$segxend=(df3$radius+0.05)*sin(df3$mid)+df3$x
      df3$segyend=(df3$radius+0.05)*cos(df3$mid)+df3$y
      if (labelposition==2)
        df3$radius=(r1+r2)/2
      df3$labelx=(df3$radius)*sin(df3$mid)+df3$x
      df3$labely=(df3$radius)*cos(df3$mid)+df3$y
    }
    else {
      df3$radius=(r1+r2)/2
      if (explodeDonut)
        df3$radius[df3$focus!=0]=df3$radius[df3$focus!=0]+df3$focus[df3$focus!=0]
      df3$labelx=df3$radius*sin(df3$mid)+df3$x
      df3$labely=df3$radius*cos(df3$mid)+df3$y
    }
    df3$segx[df3$ratio1==0]=0
    df3$segxend[df3$ratio1==0]=0
    df3$segy[df3$ratio1==0]=0
    df3$segyend[df3$ratio1==0]=0
    if (labelposition==0) {
      df3$segx[df3$ratio1<showRatioThreshold]=0
      df3$segxend[df3$ratio1<showRatioThreshold]=0
      df3$segy[df3$ratio1<showRatioThreshold]=0
      df3$segyend[df3$ratio1<showRatioThreshold]=0
    }
    df3
    del=which(df3$Freq==0)
    del
    if (length(del)>0)
      subColor<-subColor[-del]
    subColor
  }
  p<-ggplot()+theme_no_axes()+coord_fixed()
  if (is.null(maxx)) {
    r3=r2+0.3
  }
  else {
    r3=maxx
  }
  p1<-p+geom_arc_bar(aes_string(x0="x",y0="y",r0=as.character(r0),
                                r=as.character(r1),start="start1",end="end1",
                                fill=pies),alpha=pieAlpha,color=color,data=df)+
    transparent()+scale_fill_manual(values=mainCol)+
    xlim(r3*c(-1,1))+ylim(r3*c(-1,1))+guides(fill=FALSE)
  if ((labelposition==1)&(is.null(donuts))) {
    p1<-p1+geom_segment(aes_string(x="segx",y="segy",
                                   xend="segxend",yend="segyend"),data=df)+
      geom_text(aes_string(x="segxend",y="segyend",
                           label="label",hjust="hjust",vjust="vjust"),
                size=pieLabelSize,data=df,family=family)
  }
  else if ((labelposition==2)&(is.null(donuts))) {
    p1<-p1+geom_segment(aes_string(x="segx",y="segy",xend="segxend",yend="segyend"),
                        data=df[df$ratio<labelpositionThreshold,])+
      geom_text(aes_string(x="segxend",y="segyend",label="label",
                           hjust="hjust",vjust="vjust"),
                size=pieLabelSize,data=df[df$ratio<labelpositionThreshold,],
                family=family)+
      geom_text(aes_string(x="labelx",y="labely",label="label"),
                size=pieLabelSize,data=df[df$ratio>=labelpositionThreshold,],
                family=family)
  }
  else {
    p1<-p1+geom_text(aes_string(x="labelx",y="labely",label="label"),
                     size=pieLabelSize,data=df,family=family)
  }
  if (showPieName)
    p1<-p1+annotate("text",x=0,y=0,label=pies,
                    size=titlesize,family=family)
  p1<-p1+theme(text=element_text(family=family))
  if (!is.null(donuts)) {
    if (explodeDonut) {
      p3<-p+geom_arc_bar(aes_string(x0="x",y0="y",
                                    r0=as.character(r1),
                                    r=as.character(r2),
                                    start="start1",end="end1",
                                    fill="no",explode="focus"),
                         alpha=donutAlpha,color=color,data=df3)
    }
    else {
      p3<-p+geom_arc_bar(aes_string(x0="x",y0="y",
                                    r0=as.character(r1),r=as.character(r2),
                                    start="start1",end="end1",fill="no"),
                         alpha=donutAlpha,color=color,data=df3)
    }
    p3<-p3+transparent()+scale_fill_manual(values=subColor)+
      xlim(r3*c(-1,1))+ylim(r3*c(-1,1))+guides(fill=FALSE)
    p3
    if (labelposition==1) {
      p3<-p3+geom_segment(aes_string(x="segx",y="segy",
                                     xend="segxend",yend="segyend"),data=df3)+
        geom_text(aes_string(x="segxend",y="segyend",
                             label="label",hjust="hjust",vjust="vjust"),
                  size=donutLabelSize,data=df3,family=family)
    }
    else if (labelposition==0) {
      p3<-p3+geom_text(aes_string(x="labelx",y="labely",
                                  label="label"),size=donutLabelSize,data=df3,
                       family=family)
    }
    else {
      p3<-p3+geom_segment(aes_string(x="segx",y="segy",xend="segxend",yend="segyend"),
                          data=df3[df3$ratio1<labelpositionThreshold,])+
        geom_text(aes_string(x="segxend",y="segyend",label="label",
                             hjust="hjust",vjust="vjust"),
                  size=donutLabelSize,data=df3[df3$ratio1<labelpositionThreshold,],
                  family=family)+
        geom_text(aes_string(x="labelx",y="labely",label="label"),
                  size=donutLabelSize,data=df3[df3$ratio1>=labelpositionThreshold,],
                  family=family)
    }
    if (!is.null(title))
      p3<-p3+annotate("text",x=0,y=r3,label=title,
                      size=titlesize,family=family)
    else if (showDonutName)
      p3<-p3+annotate("text",x=(-1)*r3,y=r3,
                      label=donuts,hjust=0,size=titlesize,
                      family=family)
    p3<-p3+theme(text=element_text(family=family))
    grid.newpage()
    print(p1,vp=viewport(height=1,width=1))
    print(p3,vp=viewport(height=1,width=1))
  }
  else{
    p1
  }
}


