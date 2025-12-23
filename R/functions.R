#' Plot model convergence indexes
#'
#' @param Hm an Hmsc model object
#' @param beta if TRUE, plots the beta parameters (environmental covariates).
#' @param V if TRUE, plots the V parameters.
#' @param gamma if TRUE, plots the gamma parameters (traits).
#' @param omega if TRUE, plots the omega parameters (spp associations).
#' @param title character string.
#' @export
gghm_convergence <- function(Hm,
                               beta = TRUE,
                               V=FALSE,
                               gamma = FALSE,
                               omega=FALSE,
                               title = "Model Convergence"){
  requireNamespace("Hmsc")
  requireNamespace("coda")
  requireNamespace("dplyr")
  requireNamespace("tibble")

  mpost <- Hmsc::convertToCodaObject(Hm)

  d <-
    dplyr::bind_rows(
      coda::effectiveSize(mpost$Beta) |>
        tibble::as_tibble() |>
        dplyr::mutate(fit_statistic = "ess", variable = "beta"),
      coda::gelman.diag(mpost$Beta, multivariate=FALSE)$psrf|>
        tibble::as_tibble() |> dplyr::rename(value = `Point est.`) |>
        dplyr::mutate(variable = "beta", fit_statistic = "psrf")
    )

  if(V) {
    d <- d |>
      dplyr::bind_rows(
        coda::effectiveSize(mpost$V) |>
          tibble::as_tibble() |>
          dplyr::mutate(fit_statistic = "ess", variable = "V")) |>
      dplyr::bind_rows(coda::gelman.diag(mpost$V, multivariate=FALSE)$psrf|>
                  tibble::as_tibble() |> dplyr::rename(value = `Point est.`) |>
                  dplyr::mutate(variable = "V", fit_statistic = "psrf")
      )
  }

  if(gamma) {
    d <- d |>
      dplyr::bind_rows(
        coda::effectiveSize(mpost$Gamma) |>
          tibble::as_tibble() |> dplyr::mutate(fit_statistic = "ess", variable = "gamma")) |>
      dplyr::bind_rows(coda::gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf|>
                  tibble::as_tibble() |> dplyr::rename(value = `Point est.`) |>
                 dplyr::mutate(variable = "gamma", fit_statistic = "psrf")
      )
  }

  if(omega){
    sppairs = matrix(sample(x = 1:Hm$ns^2, size = 100))
    tmp = mpost$Omega[[1]]
    for (chain in 1:length(tmp)){
      tmp[[chain]] = tmp[[chain]][,sppairs]
    }

    d <- d |>
      dplyr::bind_rows(
        coda::effectiveSize(tmp) |>
          tibble::as_tibble() |> dplyr::mutate(fit_statistic = "ess", variable = "omega")) |>
      dplyr::bind_rows(coda::gelman.diag(tmp, multivariate=FALSE)$psrf|>
                  tibble::as_tibble() |> dplyr::rename(value = `Point est.`) |>
                  dplyr::mutate(variable = "omega", fit_statistic = "psrf")
      )
  }

  vline_df <- data.frame(fit_statistic = c("ess", "psrf"),
                         xintercept = c(length(mpost$Beta)*nrow(mpost$Beta[[1]]),
                                        1.01))
  if(!beta) d <- dplyr::filter(d, variable != "beta")

  d <- d |>
    dplyr::mutate(variable = stringr::str_replace_all(variable, "beta", "Environmental Covariates") |>
                    stringr::str_replace_all('gamma', "Traits") |>
                    stringr::str_replace_all('omega', "Residual Correlations"))

  ggplot2::ggplot(d, ggplot2::aes(x=value)) +
    ggplot2::geom_histogram(bins=70) +
    ggplot2::geom_vline(data = vline_df, ggplot2::aes(xintercept = xintercept), color="red", lty=2)+
    ggplot2::facet_grid(variable~fit_statistic, scales='free') +
    ggplot2::ggtitle(title)

}
#' Trace plots
#'
#' @param Hm Hmsc object
#' @param which Can be "beta", "gamma", or "v"
#' @export
gghm_traceplot <- function(Hm,
                        which = "beta"){
  requireNamespace("Hmsc")
  requireNamespace("ggmcmc")
  co <- Hmsc::convertToCodaObject(Hm)
  if(which == "beta") return(plot(co$Beta))
  if(which == "gamma") return(plot(co$Gamma))
  if(which == "v") return(plot(co$V))
}

#' Plot effective sample size
#'
#' @param Hm an Hmsc model object
#' @param beta if TRUE, plots the beta parameters (environmental covariates).
#' @param V if TRUE, plots the V parameters.
#' @param gamma if TRUE, plots the gamma parameters (traits).
#' @param omega if TRUE, plots the omega parameters (spp associations).
#' @param title character string.
#' @export
gghm_ess <- function(Hm,
                       beta = TRUE,
                       V=FALSE,
                       gamma = FALSE,
                       omega=FALSE,
                       title = "Model Convergence"){

  mpost <- Hmsc::convertToCodaObject(Hm)

  d <- coda::effectiveSize(mpost$Beta) |>
    tibble::as_tibble() |>
    dplyr::mutate(fit_statistic = "ess", variable = "beta")

  if(V) {
    d <- d |>
      dplyr::bind_rows(
        coda::effectiveSize(mpost$V) |>
          tibble::as_tibble() |>
          dplyr::mutate(fit_statistic = "ess", variable = "V"))

  }

  if(gamma) {
    d <- d |>
      dplyr::bind_rows(
        coda::effectiveSize(mpost$Gamma) |>
          tibble::as_tibble() |>
          dplyr::mutate(fit_statistic = "ess", variable = "gamma"))

  }

  if(omega){
    sppairs = matrix(sample(x = 1:Hm$ns^2, size = 100))
    tmp = mpost$Omega[[1]]
    for (chain in 1:length(tmp)){
      tmp[[chain]] = tmp[[chain]][,sppairs]
    }

    d <- d |>
      dplyr::bind_rows(
        coda::effectiveSize(tmp) |>
          tibble::as_tibble() |>
          dplyr::mutate(fit_statistic = "ess", variable = "omega"))
  }

  vline_df <- data.frame(fit_statistic = "ess",
                         xintercept = length(mpost$Beta)*nrow(mpost$Beta[[1]]))


  ggplot2::ggplot(d,ggplot2::aes(x=value)) +
    ggplot2::geom_histogram(bins=70) +
    ggplot2::geom_vline(data = vline_df,
                        ggplot2::aes(xintercept = xintercept), color="red", lty=2)+
    ggplot2::facet_grid(variable~fit_statistic, scales='free') +
    ggplot2::ggtitle(title)

}

#' Plot variance partitioning
#' @param Hm an Hmsc model object
#' @param title the title of the plot
#' @param cols a vector of colors
#' @param lut_varnames a named vector (look up table) for changing the names of variables. Useful for making things more presentable for publication.
#' @param lut_sppnames a named vector (look up table) for changing the names of species
#'
#' @examples
#' data('Hm')
#' gghm_vp(Hm)
#'
#' @export
gghm_vp <- function(Hm,
                      title = "Variance Explained",
                      cols = NULL,
                      lut_varnames = NULL,
                      lut_sppnames = NULL){
  requireNamespace("ggtext")
  requireNamespace("Hmsc")
  requireNamespace("dplyr")
  requireNamespace("tibble")
  requireNamespace("tidyr")
  VP <- Hmsc::computeVariancePartitioning(Hm)
  mpost <- Hmsc::convertToCodaObject(Hm)

  prevalence <- colSums(Hm$Y) |>
    tibble::as_tibble(rownames = "Species") |>
    dplyr::rename(prevalence = value) |>
    dplyr::arrange(dplyr::desc(prevalence))

  mf_df <- data.frame(Species = colnames(Hm$Y)) |>
    dplyr::left_join(prevalence)

  vp_df <- VP$vals|>
    tibble::as_tibble(rownames = "variable") |>
    tidyr::pivot_longer(cols = names(as.data.frame(VP$vals)),
                 names_to = "Species",
                 values_to = "value") |>
    dplyr::left_join(prevalence) |>
    stats::na.omit()

  if(is.vector(lut_varnames)) vp_df <- vp_df |> dplyr::mutate(variable = lut_varnames[variable])
  if(is.vector(lut_sppnames)) vp_df <- vp_df |> dplyr::mutate(Speices = lut_varnames[Species])

  vp_summary <- vp_df |>
    dplyr::group_by(variable) |>
    dplyr::summarise(value_pct = mean(value) * 100) |>
    dplyr::ungroup() |>
    dplyr::mutate(variable_pct = paste0(variable, " (", round(value_pct,1), "%)")) |>
    dplyr::select(-value_pct)
  vspp <- vp_df |>
    dplyr::filter(variable == dplyr::first(vp_df$variable |> unique())) |>
    dplyr::arrange(prevalence) |>
    dplyr::pull(Species)
  vp_order <- vp_df |>
    dplyr::filter(variable == dplyr::first(vp_df$variable |> unique())) |>
    dplyr::arrange(prevalence) |>
    dplyr::mutate(Species_f = factor(Species, levels = vspp)) |>
    dplyr::select(Species, Species_f)

  p_vp <- dplyr::left_join(vp_df, vp_order) |>
    dplyr::left_join(vp_summary) |>
    dplyr::mutate(variable = factor(variable),
           value = value) |>
    ggplot2::ggplot(ggplot2::aes(x=value,y=Species_f, fill = variable_pct)) +
    ggplot2::geom_bar(stat="identity", color = "black")+
    ggplot2::theme_classic() +
    ggplot2::scale_fill_discrete(name = "Variable\n (Avg Variance Explained)") +
    ggplot2::ylab("Species") +
    ggplot2::xlab("Proportion of Variance Explained") +
    ggplot2::theme(legend.position = "right",
          legend.text = ggtext::element_markdown(),
          legend.title = ggplot2::element_blank(),
          legend.justification = c(1,0),
          legend.background = ggplot2::element_rect(color="black"))

  if(is.vector(cols)) p_vp <- p_vp + ggplot2::scale_fill_manual(values = cols)
  if(is.vector(title)) p_vp <- p_vp + ggplot2::ggtitle(title)

  return(p_vp)

}

#' Plot beta posterior estimates as colored boxes
#'
#' @param Hm an Hmsc model object
#' @param order_x order the variables on the x axis
#' @param grouping_var_y group the species on the Y axis by traits
#' @param grouping_var_y_2 second grouping variable
#' @param n_gvy number of y axis grouping variables
#' @param spp_exclude a vector of species names to exclude from the figure
#' @param support_level the support level at which to make parameters visible (default = 0.89)
#' @param lut_varnames a named vector (look up table) for changing the names of environmental covariates
#' @param lut_sppnames a named vector (look up table) for changing the names of species
#' @param no_intercept exclude the intercept parameters from the plot?
#' @param title the plot title
#'
#' @export
gghm_beta <- function(Hm,
                        order_x = NA,
                        grouping_var_y = NA,
                        grouping_var_y_2 = NA,
                        n_gvy = 1,
                        spp_exclude = NA,
                        support_level = 0.89,
                        lut_varnames = NULL,
                        lut_sppnames = NULL,
                        no_intercept = TRUE,
                        title = NA){
  requireNamespace("Hmsc")
  requireNamespace("tibble")
  requireNamespace('dplyr')
  requireNamespace('ggplot2')
  requireNamespace('ggtext')
  postBeta <- Hmsc::getPostEstimate(Hm, parName = "Beta")

  covNamesNumbers <- c(TRUE, FALSE)
  covNames = character(Hm$nc)
  for (i in 1:Hm$nc) {
    sep = ""
    if (covNamesNumbers[1]) {
      covNames[i] = paste(covNames[i], Hm$covNames[i], sep = sep)
      sep = " "
    }
    if (covNamesNumbers[2]) {
      covNames[i] = paste(covNames[i], sprintf("(C%d)", i), sep = sep)
    }
  }

  means <- postBeta$mean |>
    tibble::as_tibble() |>
    tibble::rowid_to_column("env_var") |>
    dplyr::mutate(env_var = c(covNames)) |>
    tidyr::pivot_longer(cols=colnames(postBeta$mean),
                        names_to = "species", values_to = "Mean")

  supported <- postBeta$support |>
    tibble::as_tibble() |>
    tibble::rowid_to_column("env_var") |>
    dplyr::mutate(env_var = covNames) |>
    tidyr::pivot_longer(cols=colnames(postBeta$support),
                 names_to = "species",
                 values_to = "support") |>
    dplyr::filter(support > support_level | support < (1-support_level),
           env_var != "(Intercept)") |>
    dplyr::left_join(means, by = c("env_var", "species"))|>
    dplyr::mutate(sign = ifelse(Mean>0, "+", "-"))

  # calculating prevalence if it's abundance, need to reclassify
  Ypa <- Hm$Y
  Ypa[Ypa > 0] <- 1
  sp_sorted <- Ypa |>
    colSums() |>
    tibble::as_tibble(rownames = "species") |>
    dplyr::rename(prevalence = value) |>
    dplyr::arrange(prevalence) |>
    dplyr::pull(species)

  # grouping the y axis by a variable
  if(!is.na(grouping_var_y[1])){
    if(n_gvy == 1){
    sp_sorted <- Hm$TrData |>
      tibble::as_tibble(rownames = "species") |>
      dplyr::arrange(grouping_var_y) |>
      dplyr::pull(species)
    }else{
    sp_sorted <- Hm$TrData |>
      tibble::as_tibble(rownames = "species") |>
      dplyr::arrange(grouping_var_y, grouping_var_y_2) |>
      dplyr::pull(species)
      }
  }

  vp_order <-   Ypa |>
    colSums() |>
    tibble::as_tibble(rownames = "species") |>
    dplyr::rename(prevalence = value) |>
    dplyr::arrange(prevalence) |>
    dplyr::mutate(species_f = factor(species, levels = sp_sorted)) |>
    dplyr::filter(species %in% supported$species)


  supported <- supported |>
    dplyr::left_join(vp_order)#

  if(!is.na(order_x[1])){
    var_order <- tibble::tibble(var1 = order_x) |>
      dplyr::mutate(order = 1:length(order_x)) |>
      dplyr::mutate(order_f = factor(var1, levels = order_x)) |>
      dplyr::rename(env_var = var1)

    supported <- supported |>
      dplyr::left_join(var_order) |>
      dplyr::select(-env_var) |>
      dplyr::rename(env_var = order_f)

  }

  if(is.vector(lut_varnames)) supported <- supported |> dplyr::mutate(env_var = lut_varnames[env_var])
  if(is.vector(lut_sppnames)) supported <- supported |> dplyr::mutate(species = lut_varnames[species])
  if(no_intercept) supported <- supported |> dplyr::filter(env_var != "(Intercept)")
  if(!all(is.na(spp_exclude))) supported <- dplyr::filter(supported, !species %in% spp_exclude)

  p_beta <- supported |>
    ggplot2::ggplot(ggplot2::aes(x=env_var,y=stats::reorder(species_f,species))) +
    ggplot2::geom_tile(lwd=.5,ggplot2::aes(fill = Mean, color = sign)) +
    ggplot2::theme_classic()+
    ggplot2::scale_fill_steps2() +
    ggplot2::scale_color_manual(values = c(("red"), ("blue"))) +
    ggplot2::guides(color = "none")+
    ggplot2::scale_x_discrete(expand = c(0,1)) +
    ggplot2::theme(axis.text.x = ggtext::element_markdown(angle=45, vjust=1,hjust = 1),
          legend.position = "right",
          panel.grid.major.y = ggplot2::element_line(color = "grey", linetype=3),
          plot.background = ggplot2::element_rect(color="black"),
          plot.title = ggplot2::element_text(hjust = 1, face = "bold")) +
    ggplot2::xlab("Environmental Filters")+
    ggplot2::ylab("Species")

  if(!is.na(title)) p_beta <- p_beta + ggplot2::ggtitle(title)

  return(p_beta)
}

#' Plot beta estimates using PDFs
#' @param Hm an Hmsc model object
#' @param order_var order the variables on the y axis with this trait variable
#' @param grouping_var_y group the species on the Y axis by traits
#' @param grouping_var_y_2 second grouping variable
#' @param lut_gensp a named vector of species names to change them
#' @param excluded_spp a vector of species names to exclude from the figure
#' @param included_variables a vector of environmental variable names to indicate which variables to include
#' @param top_x_species a numeric, indicating to only plot the top X most prevalent species
#' @param lut_ivars a named vector of environmental variable names to change them.
#' @param group_by_trait UNDER CONSTRUCTION: a column name from the trait data (Hm$TrData) with which to sort the y axis
#' @examples
#' data("Hm")
#' gghm_beta2(Hm)
#'
#' @export
gghm_beta2 <- function(Hm,
                       order_var = 'prevalence',
                       grouping_var_y = NA,
                       grouping_var_y_2 = NA,
                       lut_gensp=NA,
                       excluded_spp = NA,
                       included_variables = NA,
                       top_x_species = NA,
                       lut_ivars = NA,
                       group_by_trait = NA){
  # todo: maybe have an option for a prevalence cutoff, or top x most prevalent
  requireNamespace('dplyr')
  requireNamespace('tidyr')
  requireNamespace('tibble')
  requireNamespace('ggplot2')
  requireNamespace('bayestestR')
  requireNamespace('ggnewscale')
  requireNamespace('ggmcmc')
  requireNamespace("stringr")

  cc<-Hmsc::convertToCodaObject(Hm)
  mbc0 <- ggmcmc::ggs(cc$Beta) |>
    tidyr::separate(col = "Parameter",
             into = c("var", 'sp'),
             sep = ", ") |>
    dplyr::mutate(var = stringr::str_remove_all(var, "B\\["),
                  var = stringr::str_remove_all(var, " \\(C\\d+\\)"),
                  sp = stringr::str_remove_all(sp, " \\(S\\d+\\)\\]"))|>
    dplyr::filter(var != "(Intercept)") |>
    dplyr::group_by(var, sp, Chain) |>
    dplyr::mutate(value = scale(value,center = F),
           sign = ifelse(value>0, "positive", "negative"),
           median_value = stats::median(value)) |>
    dplyr::filter(value<4 & value>-4) |>
    dplyr::ungroup()

  result <-  list()
  cc <- 1
  for(i in unique(mbc0$var)){
    for(j in unique(mbc0$sp)){
      result[[cc]] <- mbc0 |>
        dplyr::filter(var == i, sp == j) |>
        dplyr::pull(value) |>
        bayestestR::p_direction() |>
        dplyr::mutate(var = i, sp=j) |>
        dplyr::select(-Parameter)
      cc <- cc + 1
    }}
  pd <- result |>
    dplyr::bind_rows() |>
    tibble::as_tibble() |>
    dplyr::mutate(p_equiv = dplyr::case_when(
      pd <= 0.95 ~ 'ns',
      pd > 0.95 & pd <= 0.975 ~ '0.1',
      pd > 0.975 & pd <= 0.995 ~ '0.05',
      pd > 0.995 & pd <= 0.9995 ~ '0.01',
      pd > 0.9995 ~ '0.001'
    ))



  if(any(!is.na(included_variables))){
    mbc0 <- dplyr::filter(mbc0, var %in% included_variables)
  }

  prevalence <- Hm$Y |>
    tibble::as_tibble(rownames = "plot")
  prevalence <- prevalence |>
    tidyr::pivot_longer(cols = names(prevalence)[2:ncol(prevalence)],
                 names_to = "sp") |>
    dplyr::mutate(value = ifelse(value > 0, 1, 0)) |>
    dplyr::group_by(sp) |>
    dplyr::summarise(prevalence = sum(value),
              prev_pct = sum(value)/dplyr::n()*100) |>
    dplyr::mutate(prev_pct = ifelse(prev_pct<1, round(prev_pct,1), round(prev_pct))) |>
    dplyr::ungroup()

  if(any(!is.na(excluded_spp))){
    mbc0 <- dplyr::filter(mbc0, !sp %in% excluded_spp)
    prevalence <- dplyr::filter(prevalence, !sp %in% excluded_spp)

  }
  if(any(!is.na(lut_gensp))){
    prevalence <- dplyr::mutate(prevalence, sp = lut_gensp[sp])
  }

  mbc <- mbc0 |>
    # dplyr::left_join(Hm$TrData |>
    #                    tibble::as_tibble(rownames = "sp")) |>
    dplyr::left_join(pd)

  if(any(!is.na(lut_gensp))){
    mbc <- dplyr::mutate(mbc, sp = lut_gensp[sp])
  }

  mbc <- mbc |>
    dplyr::left_join(prevalence) |>
    dplyr::mutate(sp = paste0(sp, " (", prev_pct, ")")) |>
    dplyr::mutate(sp = stringr::str_replace_all(sp,"0.5", ".5"))

  # redo from here with the group by trait thing
  vp_order <- mbc |>
    dplyr::left_join(prevalence) |>
    dplyr::filter(var == mbc$var[1],
           Iteration ==1, Chain==1) |>
    dplyr::arrange(prevalence) # this is where the order by functional group comes in
  vp_order <- vp_order |>
    dplyr::mutate(sp_f = factor(sp, levels = vp_order$sp)) |>
    dplyr::select(sp, sp_f)

  # maybe have a whole section here where we just do this
  if(!any(is.na(group_by_trait))){
    dd <- Hm$TrData |>
      tibble::as_tibble(rownames = 'sp') |>
      dplyr::filter(!sp %in% excluded_spp) |>
      dplyr::pull(dplyr::any_of(group_by_trait))

    hline_position <- table(dd)
   }

  dp <- mbc |> dplyr::left_join(vp_order)

  if(any(!is.na(lut_ivars))){
    dp <- dplyr::mutate(dp, var = lut_ivars[var])
  }

  if(!is.na(top_x_species)){
    top_spp <- dp |>
      dplyr::select(sp, prevalence) |>
      unique() |>
      dplyr::mutate(sp_rank = length(unique(dp$sp)) -rank(prevalence)) |>
      dplyr::filter(sp_rank < top_x_species) |>
      dplyr::pull(sp)

    dp <- dplyr::filter(dp, sp %in% top_spp)
  }

  pcols <- c("ns" = "white",
             '0.1' = 'grey40',
             '0.05' = 'black',
             '0.01' = 'black',
             '0.001' = 'black')

  p <- ggplot2::ggplot(dp,
             ggplot2::aes(x=value, y = sp_f,
                  group=as.factor(Chain))) +
    ggdist::stat_slab(height=2,
                      lwd = .75,
                      ggplot2::aes(fill = ggplot2::after_stat(x>0),
                                   color = p_equiv,
                                   alpha = exp(abs(median_value))))+
    ggplot2::facet_wrap(~var,
                        scales = "free_x",
                        nrow=1,
                        ncol=length(unique(mbc$var))) +
    ggplot2::scale_alpha_continuous(range = c(0,2*(1/length(unique(mbc$Chain)))))+
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = pcols) +
    ggplot2::guides(fill="none", alpha="none", color = "none")+
    ggplot2::geom_vline(xintercept=0, col="black", lty=2) +
    ggnewscale::new_scale_fill() +
    ggplot2::xlab("Scaled Effect on Occurrence Probability") +
    ggplot2::ylab("Species or Species Group (% Prevalence)") +
    ggplot2::theme(panel.spacing.x = ggplot2::unit(-1, "lines"),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(size=12))
  if(!any(is.na(group_by_trait))) p<- p + ggplot2::geom_hline(yintercept = hline_position + 0.5)

  return(p)
}

#' Do tests on parameters
#'
#' @param Hm Hmsc object
#' @param ci Credible interval
#'
#' @export
gghm_beta_tests <- function(Hm, ci = 0.95){
  requireNamespace('dplyr')
  requireNamespace('tidyr')
  requireNamespace('tibble')
  requireNamespace('ggmcmc')
  requireNamespace('bayestestR')
  requireNamespace("stringr")

  cc<-Hmsc::convertToCodaObject(Hm)
  mbc0 <- ggmcmc::ggs(cc$Beta) |>
    tidyr::separate(col = "Parameter",
                    into = c("var", 'sp'),
                    sep = ", ") |>
    dplyr::mutate(var = stringr::str_remove_all(var, "B\\["),
                  var = stringr::str_remove_all(var, " \\(C\\d+\\)"),
                  sp = stringr::str_remove_all(sp, " \\(S\\d+\\)\\]"))|>
    dplyr::filter(var != "(Intercept)") |>
    dplyr::group_by(var, sp, Chain) |>
    dplyr::mutate(value = scale(value,center = F),
                  sign = ifelse(value>0, "positive", "negative"),
                  median_value = stats::median(value)) |>
    dplyr::filter(value<4 & value>-4) |>
    dplyr::ungroup()


    result <-  list()
    cc <- 1
    for(i in unique(mbc0$var)){
      for(j in unique(mbc0$sp)){
        result[[cc]] <- mbc0 |>
          dplyr::filter(var == i, sp == j) |>
          dplyr::pull(value) |>
          bayestestR::describe_posterior(ci = ci) |>
          dplyr::mutate(var = i, sp=j)
        result[[cc]][2,2] <- result[[cc]][1,2]
        result[[cc]] <- dplyr::filter(result[[cc]], Parameter == "Posterior") |>
          dplyr::select(-Parameter)
        cc <- cc + 1
      }}
    return(result |> dplyr::bind_rows() |> tibble::as_tibble())

}

#' plot trait covariate relationships
#'
#' @param Hm an Hmsc object
#' @param support_level make the boxes dissappear below this number
#' @param no_intercept exclude the intercept?
#' @param title plot title
#' @export
gghm_gamma <- function(Hm,
                         support_level = 0.89,
                         no_intercept = TRUE,
                         title = "Effects on Traits"){
  requireNamespace('stringr')
  requireNamespace('Hmsc')
  requireNamespace('dplyr')
  requireNamespace('tibble')
  requireNamespace('ggplot2')
  covNamesNumbers <- c(TRUE, FALSE)
  covNames = character(Hm$nc)
  for (i in 1:Hm$nc) {
    sep = ""
    if (covNamesNumbers[1]) {
      covNames[i] = paste(covNames[i], Hm$covNames[i], sep = sep)
      sep = " "
    }
    if (covNamesNumbers[2]) {
      covNames[i] = paste(covNames[i], sprintf("(C%d)", i), sep = sep)
    }
  }

  trNames = character(Hm$nt)
  trNamesNumbers = c(T,F)
  for (i in 1:Hm$nt) {
    sep = ""
    if (trNamesNumbers[1]) {
      trNames[i] = paste(trNames[i], Hm$trNames[i], sep = sep)
      sep = " "
    }
    if (trNamesNumbers[2]) {
      trNames[i] = paste(trNames[i], sprintf("(T%d)", i),
                         sep = sep)
    }
  }

  trNames <- stringr::str_remove_all(trNames, "yes") |>
    stringr::str_remove("pp")

  postGamma = Hmsc::getPostEstimate(Hm, parName="Gamma")

  means_gamma <- postGamma$mean |>
    tibble::as_tibble() |>
    tibble::rowid_to_column("env_var") |>
    dplyr::mutate(env_var = c(covNames))

  means_gammal <- means_gamma |>
    tidyr::pivot_longer(cols=names(means_gamma)[2:ncol(means_gamma)],
                        names_to = "Trait", values_to = "Mean")

  lut_gamma <- trNames
  names(lut_gamma) <- unique(means_gammal$Trait)

  supported_gamma0 <- postGamma$support |>
    tibble::as_tibble() |>
    tibble::rowid_to_column("env_var") |>
    dplyr::mutate(env_var = covNames)
  supported_gamma <- supported_gamma0 |>
    tidyr::pivot_longer(cols=names(supported_gamma0)[2:ncol(.)],
                 names_to = "Trait",
                 values_to = "Support") |>
    dplyr::filter(Support > support_level |Support< (1-support_level),
           env_var != "(Iintercept)") |>
    dplyr::left_join(means_gammal, by = c("env_var", "Trait"))|>
    dplyr::mutate(sign = ifelse(Mean>0, "+", "-"),
           Trait = lut_gamma[Trait])|>
    dplyr::filter(Trait != "(Intercept)",
                  env_var != "(Intercept)")

  p_gamma <- supported_gamma |>
    ggplot2::ggplot(ggplot2::aes(x=env_var,y=Trait, fill = Mean, color = sign)) +
    ggplot2::geom_tile(lwd=.5) +
    ggpubr::theme_pubclean()+
    ggplot2::scale_fill_steps2() +
    ggplot2::scale_color_manual(values = c(("red"), ("blue"))) +
    ggplot2::guides(color = "none")+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, vjust=1,hjust = 1),
          legend.position = "right",
          plot.background = ggplot2::element_rect(color="black"),
          plot.title = ggplot2::element_text(hjust = 1, face = "bold")) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab("Environmental Filters") +
    ggplot2::ylab("Traits")

  return(p_gamma)
}

#' plot trait covariate associations
#' @param Hm an Hmsc model object
#' @param lut_varnames a named vector used to change the variable names
#' @export
gghm_gamma2 <- function(Hm,
                          lut_varnames = NA){
  c<-Hmsc::convertToCodaObject(Hm)
  mbc <- ggmcmc::ggs(c$Gamma) |>
    tidyr::separate(.,
             col = "Parameter",
             into = c("var", "x1", "trait", "x2"),
             sep = " ") |>
    dplyr::select(-x1, -x2) |>
    dplyr::mutate(var = stringr::str_remove_all(var, "G\\["),
           trait = stringr::str_remove_all(trait, "yes"),
           trait = stringr::str_remove_all(trait, "cots"),
           trait = stringr::str_replace_all(trait, "originN", "native"),
           trait = stringr::str_to_title(trait))|>
    dplyr::filter(var != "(Intercept)", trait != "(Intercept)")

  if(!is.na(lut_varnames)){
    mbc <- dplyr::mutate(mbc, var = lut_varnames[var])
  }

  p <- ggplot2::ggplot(mbc,
              ggplot2::aes(x=value, y = trait,
                  fill=as.factor(Chain))) +
    ggdist::stat_dist_interval(alpha=0.5) +
    ggplot2::facet_wrap(~var, scales = "free_x", nrow=2,
               ncol=ceiling(length(unique(mbc$var))/2)) +
    ggplot2::theme_classic() +
    ggplot2::guides(fill="none")+
    ggplot2::geom_vline(xintercept=0, col="black", lty=2) +
    ggplot2::xlab("Effect on Occurrence Probability") +
    ggplot2::ylab("Trait") +
    ggplot2::theme(strip.text = ggtext::element_markdown())
  return(p)
}

#' create a correlation plot for species associations
#'
#' @param Hm an Hmsc object
#' @param support_level make sqaures white below this level
#' @param hc.method hierarchical clustering method
#' @param hc.order logical, whether to actually order the species according to hc.method
#' @param lut_gensp a named vector for changing species names
#' @param axis_text_colors_x a vector of one color
#' @param axis_text_colors_y a vector of one color
#' @param title the plot title
#' @export
gghm_omega <- function(Hm,
                         support_level = 0.89,
                         hc.method = "single",
                         hc.order = TRUE,
                         lut_gensp = NA,
                         axis_text_colors_x = "black",
                         axis_text_colors_y = "black",
                         title = "Residual Species Associations"){
  requireNamespace("ggcorrplot")
  requireNamespace("ggplot2")

  OmegaCor = Hmsc::computeAssociations(Hm)

  hmdf_mean <- OmegaCor[[1]]$mean |>
    as.matrix()
  if(!is.na(lut_gensp[1])){
  rownames(hmdf_mean) <- lut_gensp[rownames(hmdf_mean)]
  colnames(hmdf_mean) <- lut_gensp[colnames(hmdf_mean)]
}

  pcor1<- ggcorrplot::ggcorrplot(hmdf_mean,
                                 type = "lower",
                                 hc.order = hc.order,
                                 hc.method = hc.method,
                                 colors = c("red", "grey90", "blue"),
                                 title = title,
                                 tl.srt = 90,
                                 legend.title = "Residual\nCorrelation") +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::theme(plot.background = ggplot2::element_rect(color="black"),
          legend.position = c(0,1),
          axis.text.x = ggplot2::element_text(color = axis_text_colors_x,
                                     vjust = .05,
                                     face = "italic"),
          axis.text.y = ggplot2::element_text(color = axis_text_colors_y,
                                     face = "italic"),
          legend.justification = c(0,1),
          legend.background = ggplot2::element_rect(color="black"),
          plot.title = ggplot2::element_text(hjust = 1, face = "bold"))

  return(pcor1)
}

#' get Tjur R2 for each species
#'
#' @param Hm an Hmsc model object
#' @param title the plot title
#' @param sp_names look up table to change species names to be more presentable
#'
#' @export
gghm_r2_tjur <- function(Hm, title = "Variance Explained", sp_names = 'none'){
  requireNamespace("Hmsc")
  requireNamespace("tibble")
  requireNamespace("tidyr")
  requireNamespace("ggplot2")
  requireNamespace("ggtext")
  requireNamespace("dplyr")

  mpost <- Hmsc::convertToCodaObject(Hm)
  preds <- Hmsc::computePredictedValues(Hm)
  MF <- Hmsc::evaluateModelFit(hM=Hm, predY=preds)
  df <- MF |>
    tibble::as_tibble()
  df <- tidyr::pivot_longer(df, cols = names(df))
  spp <- colnames(Hm$Y)

    means <- df |>
      dplyr::filter(name %in% c("TjurR2", "R2")) |>
      dplyr::mutate(species = spp |> stringr::str_replace_all("_", ' ')) |>
      stats::na.omit()
    if(sp_names[1] != "none") means <- dplyr::mutate(means, species = sp_names[species])
    return(
      ggplot2::ggplot(means, ggplot2::aes(x=value, y= stats::reorder(species, value))) +
             ggplot2::geom_bar(stat = "identity") +
             ggplot2::ggtitle(title) +
             ggplot2::ylab("Species") +
             ggplot2::xlab("Tjur R<sup>2</sup>") +
             ggplot2::theme(axis.title.x = ggtext::element_markdown(),
                            axis.text.y = ggplot2::element_text(face = "italic")) +
        ggplot2::theme_bw()
      )
}


