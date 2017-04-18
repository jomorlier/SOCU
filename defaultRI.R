defaultRI <- function(repairMargin=1e-2) {
  ri = list(
    RIMODE=2   # 0: OLD RI, 1: RI w/o epsilon-feasibility, 2: RI2, the recommended case
    # 3: repairChootinan
    ,eps1=1e-5  # include all constraints not eps1-feasible into the repair mechanism
    ,eps2=1e-5  # the repair strategy looks for points in parallelepiped which are eps2-feasible in every included constraint
    ,eps3=1e-5
    # random realizations which are eps2-feasible 
    ,q=3     # draw alpha_k from uniform distribution U[0,q]
    ,mmax=1000  # draw mmax random realizations
    ,OLD=FALSE  # TRUE: activate the old repairInfeasible (before 2014-09-29) 
    ,kappa=1.5  # (only OLD) if =1.0: try to step directly to the true boundary,  
    #            if >1.0: move q bit further into the feasible region
    ,repairMargin=repairMargin    # repair only solutions whose infeasibility is less 
    # than this margin 
    ,repairOnlyFresBetter=FALSE   # if repairOnlyFresBetter=TRUE, then
    # repair only iterates with fitness < so-far-best-fitness + marFres
    ,marFres=0  # only relevant if repairOnlyFresBetter==TRUE 
  )
  return(ri);
}
