## Functions --------------------------------------------------------------------
func_main <- function(df,a,b,c)
{
  # df <-  ALl_mat
  # a <- "Benzo(a)anthracene"
  # b <- "Pirellula staleyi"
  # c <- "Syllidae_HK04"
  data_df <- func_df(df,a,b,c)
  func_step2_list <- func_step2(data_df)
  if (func_step2_list$step2_re$`Pr(>|t|)` <= 0.05) {
    func_step3_list <- func_step3(data_df)
    Step3_P <- func_step3_list$step3_re %>% filter(terms %in% "med") %>% .$`Pr(>|t|)`
    Step3_beta <- func_step3_list$step3_re %>% filter(terms %in% "iv") %>% .$Estimate %>% abs()
    #step1_beta <- func_step1_list$step1_re$Estimate %>% abs()
    #Delta_beta <- step1_beta - Step3_beta
    #if (Step3_P <= 0.05 && Delta_beta > 0) {
    if (Step3_P <= 0.05) {
      func_step1_list <- func_step1(data_df)
      #step1_beta <- func_step1_list$step1_re$Estimate %>% abs()
      func_step4_df <- func_step4(func_step2_list$fit.mediator,func_step3_list$fit.dv)
      re_df <- data.frame(
        iv = a,
        med = b,
        dv = c,
        Beta_step1_iv = func_step1_list$step1_re$Estimate,
        Beta_step1_iv_p = func_step1_list$step1_re$`Pr(>|t|)`,
        Beta_step2_iv = func_step2_list$step2_re$Estimate,
        Beta_step2_iv_p = func_step2_list$step2_re$`Pr(>|t|)`,
        Beta_step3_iv = func_step3_list$step3_re %>%
          filter(terms %in% "iv") %>% .$Estimate,
        Beta_step3_iv_p = func_step3_list$step3_re %>%
          filter(terms %in% "iv") %>% .$`Pr(>|t|)`,
        Beta_step3_med = func_step3_list$step3_re %>%
          filter(terms %in% "med") %>% .$Estimate,
        Beta_step3_med_p = func_step3_list$step3_re %>%
          filter(terms %in% "med") %>% .$`Pr(>|t|)`,
        ACME = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "ACME") %>%
          .$Estimate,
        ADE = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "ADE") %>%
          .$Estimate,
        ACME_lower_95_CI = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "ACME") %>%
          .$`95% CI Lower`,
        ACME_upper_95_CI = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "ACME") %>%
          .$`95% CI Upper`,
        ACME_P_value = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "ACME") %>%
          .$`p-value`,
        ADE_P_value = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "ADE") %>%
          .$`p-value`,
        Percentage_Medi_P_value = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "Prop. Mediated") %>%
          .$`p-value`,
        Percentage_Medi = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "Prop. Mediated") %>%
          .$Estimate,
        Percentage_Medi_lower_95_CI = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "Prop. Mediated") %>%
          .$`95% CI Lower`,
        Percentage_Medi_upper_95_CI = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "Prop. Mediated") %>%
          .$`95% CI Upper`
      )
    }else
    {
      re_df <- data.frame(
        iv = a,
        med = b,
        dv = c,
        Beta_step1_iv = "",
        Beta_step1_iv_p = "",
        Beta_step2_iv = func_step2_list$step2_re$Estimate,
        Beta_step2_iv_p = func_step2_list$step2_re$`Pr(>|t|)`,
        Beta_step3_iv = func_step3_list$step3_re %>%
          filter(terms %in% "iv") %>% .$Estimate,
        Beta_step3_iv_p = func_step3_list$step3_re %>%
          filter(terms %in% "iv") %>% .$`Pr(>|t|)`,
        Beta_step3_med = func_step3_list$step3_re %>%
          filter(terms %in% "med") %>% .$Estimate,
        Beta_step3_med_p = func_step3_list$step3_re %>%
          filter(terms %in% "med") %>% .$`Pr(>|t|)`,
        ACME = "",
        ADE = "",
        ACME_lower_95_CI = "",
        ACME_upper_95_CI = "",
        ACME_P_value = "",
        ADE_P_value = "",
        Percentage_Medi_P_value = "",
        Percentage_Medi = "",
        Percentage_Medi_lower_95_CI = "",
        Percentage_Medi_upper_95_CI = ""
      )
    }
  }else{
    #func_step3_list <- func_step3(data_df)
    re_df <- data.frame(
      iv = a,
      med = b,
      dv = c,
      Beta_step1_iv = "",
      Beta_step1_iv_p = "",
      Beta_step2_iv = func_step2_list$step2_re$Estimate,
      Beta_step2_iv_p = func_step2_list$step2_re$`Pr(>|t|)`,
      Beta_step3_iv = "",
      Beta_step3_iv_p = "",
      Beta_step3_med = "",
      Beta_step3_med_p = "",
      ACME = "",
      ADE = "",
      ACME_lower_95_CI = "",
      ACME_upper_95_CI = "",
      ACME_P_value = "",
      ADE_P_value = "",
      Percentage_Medi_P_value = "",
      Percentage_Medi = "",
      Percentage_Medi_lower_95_CI = "",
      Percentage_Medi_upper_95_CI = ""
    )
  }
  return(re_df)
}

func_df <- function(df,a,b,c)
{
  # df <-  ALl_mat
  # a <- "Nitrate Nitrogen"
  # b <- "Nitrospina gracilis"
  # c <- "Polynoidae HK03"
  data_df <- df %>%
    dplyr::select(a,b,c) %>%
    `colnames<-`(c("iv","med","dv"))
  return(data_df)
}
  
func_step1 <- function(df)
{
  #Step 1: The total effect
  fit.totaleffect=lm(dv~iv,df)
  #summary(fit.totaleffect)
  step1_re <- as.data.frame(coef(summary(fit.totaleffect))) %>%
    rownames_to_column("terms") %>%
    filter(terms %in% "iv") %>%
    mutate(model = "step1")
  return_list <- list(
    fit.totaleffect = fit.totaleffect,
    step1_re = step1_re
  )
  return(return_list)
}

func_step2 <- function(df)
{
  #Step 2: The effect of the IV onto the mediator
  fit.mediator=lm(med~iv,df)
  #summary(fit.mediator)
  step2_re <- as.data.frame(coef(summary(fit.mediator))) %>%
    rownames_to_column("terms") %>%
    filter(terms %in% "iv") %>%
    mutate(model = "step2")
  return_list <- list(
    fit.mediator = fit.mediator,
    step2_re = step2_re
  )
  return(return_list)
}

func_step3 <- function(df)
{
  #Step 3: The effect of the mediator on the dependent variable
  fit.dv=lm(dv~iv+med,df)
  #summary(fit.dv)
  step3_re <- as.data.frame(coef(summary(fit.dv))) %>%
    rownames_to_column("terms") %>%
    filter(terms %in% c("iv","med")) %>%
    mutate(model = "step3")
  return_list <- list(
    fit.dv = fit.dv,
    step3_re = step3_re
  )
  return(return_list)
}

func_step4 <- function(fit.mediator,fit.dv)
{
  #Step 4: Causal Mediation Analysis
  results = mediate(fit.mediator, fit.dv, treat='iv', mediator='med', boot=T)
  #summary(results)
  summary_tb <- extract_mediation_summary(summary(results))
  return(summary_tb)
}

## Sub functions ---------------------------------------------------------------
extract_mediation_summary <- function (x) { 
  
  clp <- 100 * x$conf.level
  isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) || 
                   (inherits(x$model.y, "glm") && x$model.y$family$family == 
                      "gaussian" && x$model.y$family$link == "identity") || 
                   (inherits(x$model.y, "survreg") && x$model.y$dist == 
                      "gaussian"))
  
  printone <- !x$INT && isLinear.y
  
  if (printone) {
    
    smat <- c(x$d1, x$d1.ci, x$d1.p)
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    
    rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")
    
  } else {
    smat <- c(x$d0, x$d0.ci, x$d0.p)
    smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
    smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
    smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
    smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
    
    rownames(smat) <- c("ACME (control)", "ACME (treated)", 
                        "ADE (control)", "ADE (treated)", "Total Effect", 
                        "Prop. Mediated (control)", "Prop. Mediated (treated)", 
                        "ACME (average)", "ADE (average)", "Prop. Mediated (average)")
    
  }
  
  colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""), 
                      paste(clp, "% CI Upper", sep = ""), "p-value")
  smat
  
}

