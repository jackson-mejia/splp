#' @title Local Projections Estimator for Smooth Local Projections.
#' @description
#' Constructs IRFs for both unsmoothed and smoothed local projections disciplined to a polynomial of your choice. The smoothing estimator is based on Barnichon & Brownlees (2019, ReStat) and extended to panel data. It can handle both OLS and IV methods for panel local projections as well as a variety of fixed effects.
#'
#' @param fml Formula object with the same notation as in the package fixest. For example, to run a panel local projection ten units out with two-way fixed effects for up to ten horizons out, one would write fml = f(y,0:10) ~ x | id + time. If there is an instrumental variable write fml = f(y,0:10) ~ 1 | id + time | x ~ z. Please consult the fixest package for more details.
#' @param paneldf A panel dataframe created with the fixest function fixest::panel(). You must specify the original dataframe in the first argument and the panel identifiers in the second. Consult the fixest package for more details, especially if your time identifier is not integer-spaced.
#' @param type Takes value "splp" for SPLP or "plp" for PLP.
#' @param targetvar The right-hand side variable for which we want to smooth the impulse response function.
#' @param q The number of knots in the B-spline basis function. Defaults to q = 3.
#' @param r Specifies the polynomial order r + 1 that splp tries to discipline the IRF to. For example, a linear IRF would take r = 2, which is the default value. The function can take a vector of polynomials. For example, set r = c(2,3,4) to get a linear IRF, a quadratic IRF, and a cubic IRF.
#' @param lambda Vector of penalty parameters lambda to try. Defaults to a logarithmic grid from 10^{-3} to 10^3 with fifty gridpoints.
#' @param boots Number of bootstraps to run. Defaults to zero.
#' @param clustervar Which variable to cluster around when doing a wild cluster bootstrap. Defaults to the individual identifier.
#' @param ci Confidence interval with default value 0.95 (for a 95% confidence interval).
#' @param m Choice for the number of random vectors in stochastic trace approximation. Default value is 1000.

#' @export
#'
#' @returns Returns an splp object as a list with three entries. The first, named plp, contains IRFs for the standard panel local projections estimator. The second, named info, contains basic info about the object (such as identifying variables). The third entry contains a list of splp objects. Each is named according to its polynomial order. If you set r = c(2,3,4), then there will be three objects named order_1, order_2, and order_3. Each underlying object has the relevant IRFs as well as a plot comparing its IRF to the standard plp IRF.
splp <- function(fml, paneldf, type, q = 3, r = 2, lambda =  10^seq(-3, 3, length.out = 100), boots = 0, clustervar = NULL,targetvar = NULL, ci = 0.95, m = 1000){
  `%notin%` <- Negate(`%in%`)
  fixest::setFixest_notes(FALSE)
  if (!exists("targetvar")){
    print("Error: Enter a variable name to be targeted in the field targetvar = ")
  }

  if(targetvar %notin% names(paneldf)){
    print("Error: The target variable is not in the dataset.")
  }

  if(is.null(clustervar)){
    plp <- suppressMessages(fixest::feols(fml ,data = paneldf, data.save = TRUE, demeaned = TRUE))
  } else{
    plp <- suppressMessages(fixest::feols(fml ,data = paneldf, data.save = TRUE, demeaned = TRUE, cluster = clustervar))
  }

  id_vars <- plp[[1]]$panel.id
  if(is.null(clustervar)){
    clustervar <- id_vars[1]
  }

  if("iv_wh" %in% names(plp[[1]])){
    targetvar <- paste0("fit_",targetvar)
  }

  coefdf <- data.frame(horizon = seq(0,4,length(plp)), cil = rep(NA,length(plp)), cih = rep(NA,length(plp)), IRF = rep(NA,length(plp)))
  for(ii in 1:length(plp)){
    vv <- summary(plp[[ii]])
    clevel <- -qnorm((1-ci)/2)
    coefdf$IRF[ii] <- vv$coefficients[which(names(vv$coefficients) == targetvar)]
    coefdf$cil[ii] <- coefdf$IRF[ii]-clevel*vv$se[which(names(vv$se) == targetvar)]
    coefdf$cih[ii] <- coefdf$IRF[ii]+clevel*vv$se[which(names(vv$se) == targetvar)]
    coefdf$horizon[ii] <- ii-1
  }
  p1 <- ggplot2::ggplot(data = coefdf) +
    ggplot2::geom_line(ggplot2::aes(x = horizon, y = IRF),linewidth = 1.25) +
    ggplot2::geom_ribbon(ggplot2::aes(x = horizon, ymin = cil, ymax = cih), fill = "grey70", alpha = 0.3) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = "IRF",
                  x = "Horizon") +
    ggplot2::theme(plot.title = ggplot2::element_text(size=18)) +
    ggplot2::theme(axis.title = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0), size = 18)) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(vjust = +1.25)) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(vjust = -0.1)) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 18, color = 'black'))+
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::theme(panel.border = ggplot2::element_blank(), axis.line = ggplot2::element_line()) +
    ggplot2::theme(axis.text=ggplot2::element_text(colour="black")) + ggplot2::theme(legend.text=ggplot2::element_text(size=15)) + ggplot2::theme(legend.title=ggplot2::element_text(size=13))
  splp_obj <- list()
  splp_obj$plp$plp <- plp
  splp_obj$plp$regplot <- p1
  splp_obj$plp$regcoefdf <- coefdf
  splp_obj$info$id_vars <- id_vars
  H <- length(plp) - 1
  splp_obj$info$fml <- H
  splp_obj$info$clustervar <- clustervar
  splp_obj$info$fml <- fml
  if(type == "plp"){
    return(splp_obj)
  } else{
    splp_obj$splp$info$targetvar <- targetvar
    baselinelist <- list()
    baselinelist$reg <- list()
    baselinelist$ydat <- list()
    baselinelist$xdat <- list()
    baselinelist$wdat <- list()
    h1 <- 0
    # construct basis f
    h.range <- h1:H
    bdeg    <- q
    splp_obj$splp$info$q <- q
    knots   <- seq(-bdeg+h1,H+bdeg,1)
    basis   <- splines::spline.des(knots, h.range, bdeg + 1, 0 * h.range ,outer.ok=TRUE)$design

    HR <- H+1-h1
    dmat <- diag(HR)
    colnames(basis) <- colnames(basis, do.NULL = FALSE, prefix = "basis")

    nwmat <- ncol(plp[[1]]$X_demeaned) - 1
    if(nwmat > 0){
      wmatcols <- do.call(cbind, replicate(nwmat, dmat, simplify=FALSE))
      controls_n <- list()
      for(j in 1:nwmat){
        tmp <- rep(NA, HR)
        for(i in 1:HR){
          tmp[i] <- paste("v", j, i-1, sep = "_")
        }
        controls_n[[j]] <- tmp
      }
      colnames(wmatcols) <- unlist(controls_n)
    }

    paneldf_dt <- data.table::as.data.table(paneldf)
    length_plp <- length(plp)
    id_vars_index <- id_vars %in% clustervar

    for (j in 1:length_plp) {
      obs_plp_j <- fixest::obs(plp[[j]])
      plp_j_y_demeaned <- plp[[j]]$y_demeaned
      plp_j_X_demeaned <- plp[[j]]$X_demeaned

      if (clustervar %in% id_vars) {
        baselinelist$ydat[[j]] <- cbind(horizon = j - 1, paneldf_dt[obs_plp_j, ..id_vars], dvar = plp_j_y_demeaned)
      } else {
        gvars <- c(clustervar, id_vars)
        baselinelist$ydat[[j]] <- cbind(horizon = j - 1, paneldf_dt[obs_plp_j, ..gvars], dvar = plp_j_y_demeaned)
      }

      if (nwmat > 0) {
        tmpb <- matrix(rep(wmatcols[j, ], nrow(plp_j_X_demeaned)), ncol = ncol(wmatcols), byrow = TRUE, dimnames = list(NULL, colnames(wmatcols)))

        dmW <- plp_j_X_demeaned[, -which(colnames(plp_j_X_demeaned) == targetvar)]
        if (is.vector(dmW)) {
          dmW <- as.matrix(dmW)
        }

        length_plp_j <- length(plp)
        for (dd in 1:nwmat) {
          tmpb[, seq(dd * length_plp_j - length_plp_j + 1, dd * length_plp_j, 1)] <- dmW[, dd] * tmpb[, seq(dd * length_plp_j - length_plp_j + 1, dd * length_plp_j, 1)]
        }

        baselinelist$wdat[[j]] <- cbind(horizon = j - 1, paneldf_dt[obs_plp_j, ..id_vars], tmpb)
      }

      dmX <- plp_j_X_demeaned[, targetvar]
      baselinelist$xdat[[j]] <- cbind(horizon = j - 1, paneldf_dt[obs_plp_j, ..id_vars], matrix(rep(basis[j, ], nrow(plp_j_X_demeaned)), ncol = ncol(basis), byrow = TRUE, dimnames = list(NULL, colnames(basis))) * dmX)
    }


    # bigx is the smoothed variable matrix, bigw is the non-smoothed, and bigy is the dependent variable
    bigx <- data.table::data.table(do.call(rbind, baselinelist$xdat))

    bigw <- data.table::data.table(do.call(rbind, baselinelist$wdat))
    bigy <- data.table::data.table(do.call(rbind, baselinelist$ydat))

    # join smoothed (by basis function) and non-smoothed data
    if(nwmat > 0){
      xw <- bigx[bigw, on = c(id_vars, "horizon")]
      xw <- with(xw, xw[order(get(id_vars[1]),get(id_vars[2]), horizon)] )
    } else{
      xw <- with(bigx, bigx[order(get(id_vars[1]),get(id_vars[2]), horizon)] )
    }
    # order correctly so that horizons and demeaned data align
    bigy <- with(bigy, bigy[order(get(id_vars[1]),get(id_vars[2]), horizon)] )
    alldat <- bigy[xw, on = c(id_vars, "horizon")]
    alldat <- with(alldat, alldat[order(get(id_vars[1]),get(id_vars[2]), horizon)] )

    # select variables
    XXw <- xw[ ,c(id_vars, "horizon"):=NULL]
    if(clustervar %in% id_vars){
      Yy <- as.matrix(bigy[,c(id_vars, "horizon"):=NULL])
    } else {
      Yy <- as.matrix(bigy[,c(gvars, "horizon"):=NULL])
    }
    # Convert explanatory variables to sparse matrix
    XXw<- Matrix::Matrix(as.matrix(XXw), sparse = TRUE)
    # Prep data for penalized regression by putting into two big matrices
    XX <- Matrix::t(XXw)%*%XXw
    XY <- Matrix::t(XXw)%*%Yy

    # penalty

    XS <- ncol(basis)
    splp_obj$splp$info$gcv_vec <- rep(NA,length(r))
    for(rr in 1:length(r)){
      P <- matrix(0,ncol(XXw),ncol(XXw))
      P <- Matrix::Matrix(P,sparse=TRUE)
      D   <- diag(XS)
      for (k in seq_len(r[rr]+1)) D <- diff(D)

      P[1:XS,1:XS] <- t(D) %*% D
      # impulse response horizon
      ir    <- matrix(0,H+1,length(lambda))
      #see barnichon and brownlees for theta
      theta <- matrix(0,ncol(XXw),length(lambda))
      mul   <- matrix(0,HR,length(lambda))
      TS  <-     (length(unique(plp[[1]]$data[,id_vars[2]])))^2
      xss <- nrow(Yy)
      gcv <- rep(NA,length(lambda))

      for( i in 1:length(lambda) ){
        A         <- XX + lambda[i]*P*TS
        b         <- XY
        theta[,i] <- as.vector( Matrix::solve( A , b ) )

        beta  <- theta[1:XS,i]
        #smooth coefficients with basis function
        mul[,i]   <- as.matrix(basis) %*% as.vector(beta)
        # extract smoothed impulse rsponse
        ir[(h1+1):(H+1),i]   <- mul[,i]

        p = ncol(P)
        Z <- matrix(rnorm(p * m), nrow = p, ncol = m)
        A <- XX + lambda[i] * P  # Assuming XX and P are 13 x 13
        AZ <- A %*% Z
        # Solve the linear system A %*% invAZ = Z for invAZ
        invAZ <- solve(A, Z)
        # Estimate the trace of the hat matrix
        trace_H_approx <- sum(colSums(Z * invAZ)) / m  # Sum of element-wise products divided by m

        # Suppose y and XY (X'y) are defined for calculating residuals and MSE
        y_hat <- XXw %*% Matrix::solve(A, XY)  # Solving for ridge regression coefficients
        residuals <- Yy - y_hat
        mse <- residuals^2  # Mean squared error
        mse <- mean(mse@x)
        n <- nrow(Yy)  # Number of observations, assuming y and X are compatible
        gcv[i] <- mse /( ((1 - log(trace_H_approx )/ n )^2))

      }

      lambda_opt_idx = which(gcv == min(gcv))
      if(lambda[lambda_opt_idx] == max(lambda)){
        lambda = 10^seq(3, 6, length.out = 100)
        ir    <- matrix(0,H+1,length(lambda))
        #see barnichon and brownlees for theta
        theta <- matrix(0,ncol(XXw),length(lambda))
        mul   <- matrix(0,HR,length(lambda))
        TS  <-     (length(unique(plp[[1]]$data[,id_vars[2]])))^2
        xss <- nrow(Yy)
        gcv <- rep(NA,length(lambda))

        for( i in 1:length(lambda) ){
          A         <- XX + lambda[i]*P*TS
          b         <- XY
          theta[,i] <- as.vector( Matrix::solve( A , b ) )

          beta  <- theta[1:XS,i]
          #smooth coefficients with basis function
          mul[,i]   <- as.matrix(basis) %*% as.vector(beta)
          # extract smoothed impulse rsponse
          ir[(h1+1):(H+1),i]   <- mul[,i]

          p = ncol(P)
          Z <- matrix(rnorm(p * m), nrow = p, ncol = m)
          A <- XX + lambda[i] * P  # Assuming XX and P are 13 x 13
          AZ <- A %*% Z
          # Solve the linear system A %*% invAZ = Z for invAZ
          invAZ <- solve(A, Z)
          # Estimate the trace of the hat matrix
          trace_H_approx <- sum(colSums(Z * invAZ)) / m  # Sum of element-wise products divided by m

          # Suppose y and XY (X'y) are defined for calculating residuals and MSE
          y_hat <- XXw %*% Matrix::solve(A, XY)  # Solving for ridge regression coefficients
          residuals <- Yy - y_hat
          mse <- residuals^2  # Mean squared error
          mse <- mean(mse@x)
          n <- nrow(Yy)  # Number of observations, assuming y and X are compatible
          gcv[i] <- mse /( ((1 - log(trace_H_approx )/ n )^2))

        }
      }

      splp_obj_tmp <- splp_obj
      lambda_opt_idx = which(gcv == min(gcv))
      splp_obj_tmp$splp$info$gcv = gcv[lambda_opt_idx]
      gcvtmp <- splp_obj_tmp$splp$info$gcv
      splp_obj_tmp$splp$info$idx.opt <- lambda_opt_idx
      splp_obj_tmp$splp$info$lambda <- lambda
      splp_obj_tmp$splp$info$lambda_opt <- lambda[lambda_opt_idx]
      splp_obj_tmp$splp$ir.opt  <- cbind(horizon = seq(h1,H,1), IRF = ir[ , splp_obj_tmp$splp$info$idx.opt ])
      colnames(splp_obj_tmp$splp$ir.opt) <- c("horizon", "IRF")

      splp_obj_tmp$splp$info$alldat <- alldat
      splp_obj_tmp$splp$info$TS     <- TS
      splp_obj_tmp$splp$info$Y      <- Yy
      splp_obj_tmp$splp$info$X      <- XXw
      splp_obj_tmp$splp$info$theta  <- theta
      splp_obj_tmp$splp$info$mul    <- mul
      splp_obj_tmp$splp$info$P      <- P
      splp_obj_tmp$splp$info$ir     <- cbind(horizon = seq(h1,H,1), IRF = ir)
      splp_obj_tmp$splp$info$r <- r[rr]
      ## Create matrix of residuals with id & reg
      if(clustervar %in% id_vars){
        cnames <- c(id_vars, "horizon","dvar")
      } else{
        cnames <- c(gvars, "horizon","dvar")
      }

      ytmp <- as.matrix(Yy - XXw %*% theta[,lambda_opt_idx])
      xtmp <- as.matrix(XXw %*% theta[,lambda_opt_idx])
      residmat <- cbind(ytmp,xtmp , alldat[,..cnames])
      colnames(residmat)[1:2] <- c("residual","yhat")
      splp_obj_tmp$splp$info$residmat <- residmat
      splp_obj_tmp$splp$info$basis <- basis
      if(boots > 0){
        splp_obj_tmp <- boot_splp(splp_obj_tmp,boots)
      }
      splp_obj_tmp <- splp_ci_out(splp_obj_tmp,ci,boots)
      tmpl <- list(splp_obj_tmp$smoothplot,splp_obj_tmp$comb_plot ,splp_obj_tmp$splp$info$smoothed_coefdf,splp_obj_tmp$splp$info)
      names(tmpl) <- c("plot","comp_plot","irf","info")
      splp_obj$splp$tmpl <- tmpl
      names(splp_obj$splp)[names(splp_obj$splp) == "tmpl"] <- paste("order",r[rr],sep="_")

      splp_obj$splp$info$gcv_vec[rr] <- gcvtmp
    }

    if(length(r) > 1){
      opt_idx <- which(splp_obj$splp$info$gcv_vec == min(splp_obj$splp$info$gcv_vec))

      splp_obj$splp$info$opt_poly_idx <- opt_idx
    }

    return(splp_obj)
  }
}


#' Wild cluster bootstrap function.
#' @param obj An splp object.
#' @param boots The number of bootstraps.
#' @returns Returns an H x boots data.frame, where H is the maximum horizon length for the IRF. If the original splp object is named tmp and the polynomial order is three, then the returned bootstrap df will be located under tmp$splp$order_3$info$irboot.
boot_splp <- function(obj,boots){
  ir_boot <- list()
  basis <- obj$splp$info$basis
  plptmp <- obj$plp$plp
  id_vars <- obj$info$id_vars
  clustervar <- obj$info$clustervar
  gvars <- if(clustervar %in% id_vars){
    id_vars
  } else{
    c(clustervar,id_vars)
  }
  targetvar <- obj$splp$info$targetvar
  q <- obj$splp$info$q
  pb <- txtProgressBar(min = 0, max = boots, style = 3)

  for(bb in 1:boots){

    if("iv_wh" %in% names(plptmp[[1]])){
      fixef_vars <- plptmp[[1]]$fixef_vars
      for(ii in 1:length(plptmp)){
        endovars <- plptmp[[ii]]$iv_endo_names
        tmpmat <- matrix(NA,ncol = length(endovars),nrow = nrow(plptmp[[ii]]$X_demeaned),
                         dimnames = list(NULL,endovars))
        for(jj in 1:length(endovars)){
          fsiv <- plptmp[[ii]]$iv_first_stage[[jj]]
          tmp <- data.table::setDT(cbind(plptmp[[ii]]$data[fixest::obs(fsiv),c(unique(c(clustervar,fixef_vars,id_vars)),endovars[jj])], residual = fsiv$residuals))
          tmp <- tmp[, nresid := rnorm(1,mean = 0, sd = 1), by = clustervar]
          tmp <- tmp[, endovars[jj] := get(endovars[jj]) - residual + nresid*residual]
          tmp <- fixest::demean(X = tmp[,get(endovars[jj])], f = tmp[,..fixef_vars])
          tmpmat[,jj] <- fitted.values(lm(tmp ~ plptmp[[ii]]$iv_first_stage[[jj]]$X_demeaned))
        }
        idx <- which(colnames(plptmp[[ii]]$X_demeaned) %in% plptmp[[ii]]$iv_endo_names_fit)
        plptmp[[ii]]$X_demeaned[,idx] <- tmpmat
      }

      baselinelist <- list()
      baselinelist$reg <- list()
      baselinelist$xdat <- list()
      baselinelist$wdat <- list()
      H <- length(plptmp) - 1
      h1 <- 0
      # construct basis f
      h.range <- h1:H
      bdeg    <- q
      knots   <- seq(-bdeg+h1,H+bdeg,1)
      basis   <- splines::spline.des(knots, h.range, bdeg + 1, 0 * h.range ,outer.ok=TRUE)$design

      HR <- H+1-h1
      dmat <- diag(HR)
      colnames(basis) <- colnames(basis, do.NULL = FALSE, prefix = "basis")

      nwmat <- ncol(plptmp[[1]]$X_demeaned) - 1
      if(nwmat > 0){
        wmatcols <- do.call(cbind, replicate(nwmat, dmat, simplify=FALSE))
        controls_n <- list()
        for(j in 1:nwmat){
          tmp <- rep(NA, HR)
          for(i in 1:HR){
            tmp[i] <- paste("v", j, i-1, sep = "_")
          }
          controls_n[[j]] <- tmp
        }
        colnames(wmatcols) <- unlist(controls_n)
      }

      for(j in 1:length(plptmp)){

        if(nwmat > 0){
          tmpb =  matrix(rep(wmatcols[j,],nrow(plp[[j]]$X_demeaned)),ncol = ncol(wmatcols),byrow = TRUE,  dimnames = list(NULL, colnames(wmatcols)))

          dmW <- plp[[j]]$X_demeaned[,-which(colnames(plp[[j]]$X_demeaned) == targetvar)]
          if(is.vector(dmW)){
            dmW <- as.matrix(dmW)
          }
          for(dd in 1:nwmat){
            tmpb[,seq(dd*length(plp)-length(plp) + 1,dd*length(plp),1)] = dmW[,dd]*tmpb[,seq(dd*length(plp)-length(plp) + 1,dd*length(plp),1)]
          }
          baselinelist$wdat[[j]] <- cbind(horizon = j-1, as.data.frame(paneldf)[fixest::obs(plp[[j]]),id_vars], tmpb)
        }
        dmX <- plptmp[[j]]$X_demeaned[,targetvar]
        baselinelist$xdat[[j]] <- cbind(horizon = j-1, as.data.frame(plptmp[[j]]$data)[fixest::obs(plptmp[[j]]),id_vars], matrix(rep(basis[j,],nrow(plptmp[[j]]$X_demeaned)),ncol = ncol(basis),byrow = TRUE,  dimnames = list(NULL, colnames(basis)))*dmX)
      }
      # bigx is the smoothed variable matrix, bigw is the non-smoothed, and bigy is the dependent variable
      bigx <- data.table::data.table(do.call(rbind, baselinelist$xdat))
      bigw <- data.table::data.table(do.call(rbind, baselinelist$wdat))

      # join smoothed (by basis function) and non-smoothed data
      if(nwmat > 0){
        xw <- bigx[bigw, on = c(id_vars, "horizon")]
        xw <- with(xw, xw[order(get(id_vars[1]),get(id_vars[2]), horizon)] )
      } else{
        xw <- with(bigx, bigx[order(get(id_vars[1]),get(id_vars[2]), horizon)] )
      }
      # select variables
      XXw <- xw[ ,c(id_vars, "horizon"):=NULL]
      # Convert explanatory variables to sparse matrix
      XXw<- Matrix::Matrix(as.matrix(XXw), sparse = TRUE)
    } else{
      XXw <- obj$splp$info$X
    }
    nbigy <- data.table::setDT(obj$splp$info$residmat)
    nbigy <- nbigy[, nresid := rnorm(1,mean = 0, sd = 1), by = clustervar]
    nbigy <- nbigy[,dvarn := dvar - residual + nresid*residual]

    # order correctly so that horizons and demeaned data align
    bigy <- with(nbigy, nbigy[order(get(id_vars[1]),get(id_vars[2]), horizon)] )
    Yy <- as.matrix(bigy[,"dvarn"])
    XX <- Matrix::t(XXw)%*%XXw
    XY <- Matrix::t(XXw)%*%Yy

    r <- obj$splp$info$r

    # penalty
    P <- matrix(0,ncol(XXw),ncol(XXw))
    P <- Matrix::Matrix(P,sparse=TRUE)
    XS <- ncol(basis)
    D   <- diag(XS)
    for (k in seq_len(r+1)) D <- diff(D)

    P[1:XS,1:XS] <- t(D) %*% D
    # impulse response horizon
    lambda <- obj$splp$info$lambda_opt
    #see barnichon and brownlees for theta
    TS  <-     obj$splp$info$TS
    A         <- XX + lambda*TS*P
    b         <- XY
    theta <- as.vector( Matrix::solve( A , b ) )

    beta  <- theta[1:XS]
    #smooth coefficients with basis function
    mul   <- as.matrix(basis) %*% as.vector(beta)
    # extract smoothed impulse rsponse
    ir_boot[[bb]] <- data.frame(horizon = seq(0,(length(plptmp)-1)), ir = mul)



    # Create progress bar string
    setTxtProgressBar(pb, bb)

  }
  obj$splp$info$irboot <- rbindlist(ir_boot)
  return(obj)
}


#' Plots IRFs and computes confidence interval from bootstraps.
#' @param obj An splp object.
#' @param ci A confidence interval.
#' @param boots The number of bootstraps.
#' @returns Returns a smoothed IRF plot and a plot comparing the smoothed to the unsmoothed plots.  If the original splp object is named tmp and the polynomial order is three, then the returned combined plot will be located under tmp$splp$order_3[[2]] and the smoothed IRF alone will be under tmp$splp$order_3[[2]]$smoothplot. If boots = 0, then only point estimates will be returned.
splp_ci_out <- function(obj,ci,boots){
  smoothed_coefdf <- data.frame(obj$splp$ir.opt)
  if(boots > 0){
    alphal = (1-ci)/2
    alphah = ci + alphal
    smoothed_coefdf <- obj$splp$info$irboot[, .(cil = quantile(ir, alphal),
                                                cih = quantile(ir, alphah)), by = horizon]
    smoothed_coefdf <- merge(smoothed_coefdf,data.table(obj$splp$ir.opt))
    obj$splp$info$smoothed_coefdf <- smoothed_coefdf

    obj$smoothplot <- ggplot2::ggplot(data = smoothed_coefdf) +
      ggplot2::geom_line(ggplot2::aes(x = horizon, y = IRF),linewidth = 1.25) +
      ggplot2::geom_ribbon(ggplot2::aes(x = horizon, ymin = cil, ymax = cih), fill = "grey70", alpha = 0.3) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(y = "IRF",
                    x = "Horizon") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(size=18)) +
      ggplot2::theme(axis.title = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0), size = 18)) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(vjust = +1.25)) +
      ggplot2::theme(axis.title.x = ggplot2::element_text(vjust = -0.1)) +
      ggplot2::theme(axis.text = ggplot2::element_text(size = 18, color = 'black'))+
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::theme(panel.border = ggplot2::element_blank(), axis.line = ggplot2::element_line()) +
      ggplot2::theme(axis.text=ggplot2::element_text(colour="black")) +
      ggplot2::theme(legend.text=ggplot2::element_text(size=15)) + ggplot2::theme(legend.title=ggplot2::element_text(size=13))
    obj$comb_plot <- ggplot2::ggplot(data = obj$plp$regcoefdf) +
      ggplot2::geom_ribbon(ggplot2::aes(x = horizon, ymin = cil, ymax = cih), fill = "0072B2", alpha = 0.3) +
      ggplot2::geom_line(ggplot2::aes(x = horizon, y = IRF, color = "Unsmoothed"),linewidth = 1.25) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::geom_ribbon(ggplot2::aes(x = horizon, ymin = cil, ymax = cih), fill = "D55E00", alpha = 0.3, data = smoothed_coefdf) +
      ggplot2::geom_line(ggplot2::aes(x = horizon, y = IRF, color = "Smoothed"),linewidth = 1.25, data = smoothed_coefdf) +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_manual(values = c("Unsmoothed" = "#0072B2", "Smoothed" = "D55E00")) +
      ggplot2::labs(y = "IRF",
                    x = "Horizon",
                    color = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(size=18)) +
      ggplot2::theme(axis.title = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0), size = 18)) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(vjust = +1.25)) +
      ggplot2::theme(axis.title.x = ggplot2::element_text(vjust = -0.1)) +
      ggplot2::theme(axis.text = ggplot2::element_text(size = 18, color = 'black'))+
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::theme(panel.border = ggplot2::element_blank(), axis.line = ggplot2::element_line()) +
      ggplot2::theme(axis.text=ggplot2::element_text(colour="black")) +
      ggplot2::theme(legend.text=ggplot2::element_text(size=12),
                     legend.position = "top")
  } else{
    obj$smoothplot <- ggplot2::ggplot(data = smoothed_coefdf) +
      ggplot2::geom_line( ggplot2::aes(x = horizon, y = IRF),linewidth = 1.25) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::theme_minimal() +
      ggplot2::labs(y = "IRF",
                    x = "Horizon") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title =  ggplot2::element_text(size=18)) +
      ggplot2::theme(axis.title =  ggplot2::element_text(margin =  ggplot2::margin(t = 0, r = 10, b = 0, l = 0), size = 18)) +
      ggplot2::theme(axis.title.y =  ggplot2::element_text(vjust = +1.25)) +
      ggplot2::theme(axis.title.x =  ggplot2::element_text(vjust = -0.1)) +
      ggplot2::theme(axis.text =  ggplot2::element_text(size = 18, color = 'black'))+
      ggplot2::theme(panel.grid.major =  ggplot2::element_blank(), panel.grid.minor =  ggplot2::element_blank()) +
      ggplot2::theme(panel.border =  ggplot2::element_blank(), axis.line =  ggplot2::element_line()) +
      ggplot2::theme(axis.text= ggplot2::element_text(colour="black")) +  ggplot2::theme(legend.text= ggplot2::element_text(size=15)) +  ggplot2::theme(legend.title= ggplot2::element_text(size=13))

    obj$comb_plot <- ggplot2::ggplot(data = obj$plp$regcoefdf) +
      ggplot2::geom_line( ggplot2::aes(x = horizon, y = IRF, color = "Unsmoothed"),linewidth = 1.25) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::geom_line( ggplot2::aes(x = horizon, y = IRF, color = "Smoothed"),linewidth = 1.25, data = smoothed_coefdf) +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_manual(values = c("Unsmoothed" = "0072B2", "Smoothed" = "D55E00")) +
      ggplot2::labs(y = "IRF",
                    x = "Horizon",
                    color = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title =  ggplot2::element_text(size=18)) +
      ggplot2::theme(axis.title =  ggplot2::element_text(margin =  ggplot2::margin(t = 0, r = 10, b = 0, l = 0), size = 18)) +
      ggplot2::theme(axis.title.y =  ggplot2::element_text(vjust = +1.25)) +
      ggplot2::theme(axis.title.x =  ggplot2::element_text(vjust = -0.1)) +
      ggplot2::theme(axis.text =  ggplot2::element_text(size = 18, color = 'black'))+
      ggplot2::theme(panel.grid.major =  ggplot2::element_blank(), panel.grid.minor =  ggplot2::element_blank()) +
      ggplot2::theme(panel.border =  ggplot2::element_blank(), axis.line =  ggplot2::element_line()) +
      ggplot2::theme(axis.text= ggplot2::element_text(colour="black")) +
      ggplot2::theme(legend.text= ggplot2::element_text(size=12),
                     legend.position = "top")
  }

  return(obj)
}
