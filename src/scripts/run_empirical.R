library(here)
source(here("src/R", "config.R"))         



run_empirical <- function() {
  
  dat <- import_data(drop_first_col = TRUE, n_assets = 7)
  
  
  #cfg_variants <- list(list(label = "test"))   
  cfg_variants <- list(
    #list(label = "std", q_flip = 0.2),
    #list(label = "tailW_tauL0.2_taileps0.5",  use_weighted_ll = TRUE,
    #     tauL = 0.2, joint_k = 2L, tail_eps = 0.50, q_flip = 0.2),
    #list(label = "tailW_tauL0.1_taileps0.3",  use_weighted_ll = TRUE,
    #     tauL = 0.1, joint_k = 2L, tail_eps = 0.30, q_flip = 0.2),
    # list(label = "tip_w252_M2000_tipk10",  use_tail_informed_prior = TRUE, W = 252L, M=2000L, tip_k=25),#, q_flip = 0.2)
    # list(label = "tip_w252_M2000_tipk5",  use_tail_informed_prior = TRUE, W = 252L, M=2000L, tip_k=13),
    # list(label = "tip_w252_M2000_tipk15",  use_tail_informed_prior = TRUE, W = 252L, M=2000L, tip_k=38),
    # 
    # list(label = "tip_w126_M2000_tipk10",  use_tail_informed_prior = TRUE, W = 126L, M=2000L, tip_k=13),#, q_flip = 0.2)
    # list(label = "tip_w126_M2000_tipk5",  use_tail_informed_prior = TRUE,  W = 126L, M=2000L, tip_k=6),
    # list(label = "tip_w126_M2000_tipk15",  use_tail_informed_prior = TRUE, W = 126L, M=2000L, tip_k=19),
    # 
    # list(label = "tip_w504_M2000_tipk10",  use_tail_informed_prior = TRUE, W = 504L, M=2000L, tip_k=50),#, q_flip = 0.2)
    # list(label = "tip_w504_M2000_tipk5",  use_tail_informed_prior = TRUE, W = 504L, M=2000L, tip_k=25),
    # list(label = "tip_w504_M2000_tipk15",  use_tail_informed_prior = TRUE, W = 504L, M=2000L, tip_k=76)
    
    list(label = "tip_w252_M2000_tip",  use_tail_informed_prior = TRUE, W = 252L, M=2000L),#, q_flip = 0.2)
    list(label = "tip_w126_M2000_tip",  use_tail_informed_prior = TRUE, W = 126L, M=2000L),#, q_flip = 0.2)
    list(label = "tip_w504_M2000_tip",  use_tail_informed_prior = TRUE, W = 504L, M=2000L)#, q_flip = 0.2)

  )
  
 v=cfg_variants[[1]]
  for (v in cfg_variants) {
    cfg <- modifyList(build_cfg(ncol(dat$U), trunc_tree=4), v[ setdiff(names(v),"label") ])
    cfg$label <- v$label %||% "cfg"
    set.seed(cfg$seed)
    
    data = dat
    res <- smc_full(dat, cfg)
    res$cfg <- cfg
    res$alg <- if (isTRUE(cfg$use_weighted_ll)) "tail_weighted" else "standard"
    save_result(res, here("empirical_results"))
  }
}

run_empirical()
















