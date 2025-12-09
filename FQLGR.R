############################### funciones auxiliares para evaluar #####################################
# 
# predecir_probabilidades_validacion_en_base_fpca <- function(res, fd_valid_centered, estimationBasis) {
#   change_basis <- inprod(fd_valid_centered$basis, estimationBasis)
#   A_valid <- t(fd_valid_centered$coefs)
#   A_psi_valid <- A_valid %*% change_basis
#   
#   beta  <- as.vector(res$beta)
#   gamma <- res$gamma
#   alpha <- res$alpha
#   
#   lin_pred  <- as.vector(A_psi_valid %*% beta)
#   quad_pred <- rowSums((A_psi_valid %*% gamma) * A_psi_valid)
#   y_prob_valid <- 1 / (1 + exp(-(alpha + lin_pred + quad_pred)))
#   return(y_prob_valid)
# }

distancia_euclidea <- function(matriz1, matriz2){
  
  proj1 <- matriz1 %*% t(matriz1)
  proj2 <- matriz2 %*% t(matriz2)
  
  err_subspace <- norm(proj1 - proj2, "F")
  return(err_subspace)
}

dist_euclidea <- function(a, b) {
  sqrt(sum((a - b)^2))
}

predecir_probabilidades_validacion_en_base_fpca <- function(res, fd_valid_centered, estimationBasis) {
  scores <- inprod(fd_valid_centered, estimationBasis)
  
  beta  <- as.vector(res$beta)
  gamma <- res$gamma
  alpha <- res$alpha
  
  lin_pred  <- as.vector(scores %*% beta)
  quad_pred <- rowSums((scores %*% gamma) * scores)
  y_prob_valid <- 1 / (1 + exp(-(alpha + lin_pred + quad_pred)))
  return(y_prob_valid)
}
predecir_probabilidades_validacion <- function(res, fd_valid_centered, estimationBasis, modo = c("PCA", "NoPCA", "PLSR")) {
  modo <- match.arg(modo)
  psi <- inprod(estimationBasis, estimationBasis)
  
  A_valid <- t(fd_valid_centered$coefs)
  A_psi_valid <- A_valid %*% psi
  
  beta <- res$beta
  gamma <- res$gamma
  alpha <- res$alpha
  
  
  
  lin_pred <- as.vector(A_psi_valid %*% beta)
  quad_pred <- rowSums((A_psi_valid %*% gamma) * A_psi_valid)
  y_prob_valid <- 1 / (1 + exp(-(alpha + lin_pred + quad_pred)))
  return(y_prob_valid)
}

  evaluar_modelos_distintas_bases <- function(
    data = list(),
    basises = list(),
    densidades_grillas = c(),
    seed = 123,
    Lcoef = c(2),
    params_Upsilon0 = list(),
    params_Upsilon0_equiv = list(),
    params_Upsilon1 = list()
) {
  
  library(pROC)
  set.seed(seed)
  resultados <- data.frame()
  distancias_df <- tibble::tibble(
    nbasis = numeric(),
    modelo1 = character(),
    modelo2 = character(),
    dist_alpha = numeric(),
    dist_beta = numeric(),
    dist_gamma = numeric()
  )
  
  # Utilidad: extrae un parámetro, o usa un valor por defecto
  get_param <- function(lst, name, default) {
    if (!is.null(lst[[name]])) return(lst[[name]])
    return(default)
  }
  
  for (basis in basises) {
    cond_psi = kappa(inprod(basis, basis))
    
    cat("====================================\n")
    cat("nbasis =", basis$nbasis, "\n")
    cat("n° cond psi: ", cond_psi )
    cat("====================================\n")
    
    X_train <- data$X_train
    X_valid <- data$X_valid
    y_train <- data$y_train
    y_valid <- data$y_valid
    
    rangevals <- basis$rangeval
    
    for (densidad in densidades_grillas) {
      
      t <- seq(rangevals[1], rangevals[2], length.out = densidad)

      harmacclLfd <- vec2Lfd(Lcoef, rangevals)
      fdPar_obj <- fdPar(basis)
      
      # ---------- Suavizado ----------
      fd_train <- smooth.basis(argvals = t, y = t(X_train), fdParobj = fdPar_obj)$fd
      fd_train_centered <- center.fd(fd_train)

      fd_valid <- smooth.basis(argvals = t, y = t(X_valid), fdParobj = fdPar_obj)$fd
      fd_valid_centered <- center.fd(fd_valid) 
      # ---------- Inicialización ----------
      nb <- basis$nbasis
      beta_init  <- rnorm(nb, sd = 0.1)
      gamma_init <- matrix(rnorm(nb^2, sd = 0.1), nrow = nb)
      
      var_fd_train <- var.fd(fd_train)
      eval_var_fd_train <- eval.bifd(t,t,var_fd_train)
      cond_var_fd_train <- kappa(eval_var_fd_train)
      # ================================
      #   MODELO FPCA — Upsilon0
      # ================================
      
      modelo_fpca <- fpca_Upsilon0(
        fd_centered = fd_train_centered,
        y = y_train,
        beta = beta_init,
        gamma = gamma_init,
        alpha = 0,
        step_gradient = get_param(params_Upsilon0, "step_gradient", 0.005),
        iterations   = get_param(params_Upsilon0, "iterations", 15000),
        tol          = get_param(params_Upsilon0, "tol", 1e-8),
        var_threshold = get_param(params_Upsilon0, "var_threshold", 0.8),
        basis = basis,
        LdPenalization = harmacclLfd,
        lambda_lin  = get_param(params_Upsilon0, "lambda_lin", 1e-4),
        lambda_quad = get_param(params_Upsilon0, "lambda_quad", 0),
        modelo_quad = get_param(params_Upsilon0, "modelo_quad", TRUE),
        verbose = FALSE
      )
      
      # ================================
      #   MODELO FPCA EQUIVALENTE — Upsilon0.2
      # ================================
      
      modelo_fpca_equiv <- fpca_upsilon0_equiv(
        fd_centered = fd_train_centered,
        y = y_train,
        beta = beta_init,
        gamma = gamma_init,
        alpha = 0,
        step_gradient = get_param(params_Upsilon0, "step_gradient", 0.005),
        iterations   = get_param(params_Upsilon0, "iterations", 15000),
        tol          = get_param(params_Upsilon0, "tol", 1e-8),
        var_threshold = get_param(params_Upsilon0, "var_threshold", 0.8),
        basis = basis,
        LdPenalization = harmacclLfd,
        lambda_lin  = get_param(params_Upsilon0, "lambda_lin", 1e-4),
        lambda_quad = get_param(params_Upsilon0, "lambda_quad", 0),
        modelo_quad = get_param(params_Upsilon0, "modelo_quad", TRUE),
        verbose = FALSE
      )
      
      
      # ================================
      #   MODELO PCA COEF — Upsilon1
      # ================================
      
      res_pca_quad <- pca_coef_Upsilon1(
        fd_centered = fd_train_centered,
        y = y_train,
        beta = beta_init,
        gamma = gamma_init,
        alpha = 0,
        step_gradient = get_param(params_Upsilon1, "step_gradient", 0.005),
        iterations   = get_param(params_Upsilon1, "iterations", 10000),
        tol          = get_param(params_Upsilon1, "tol", 1e-8),
        var_threshold = get_param(params_Upsilon1, "var_threshold", 0.8),
        basis = basis,
        LdPenalization = harmacclLfd,
        lambda_lin  = get_param(params_Upsilon1, "lambda_lin", 0),
        lambda_quad = get_param(params_Upsilon1, "lambda_quad", 0),
        modelo_quad = get_param(params_Upsilon1, "modelo_quad", TRUE),
        verbose = FALSE
      )
      
      
      # =======================
      #   AUC
      # =======================
      y_prob_fpca <- predecir_probabilidades_validacion_en_base_fpca(
        modelo_fpca, fd_valid_centered, modelo_fpca$fpca_model$harmonics)
      
      y_prob_fpca_equiv <- predecir_probabilidades_validacion(
        modelo_fpca_equiv, fd_valid_centered, basis, modo = "PCA")
      
      y_prob_pca <- predecir_probabilidades_validacion(
        res_pca_quad, fd_valid_centered, basis, modo = "PCA")
      
      # =======================
      # Distancias entre parametros
      # =======================
      modelos <- list(
        fpca       = modelo_fpca,
        equiv      = modelo_fpca_equiv,
        pca_quad   = res_pca_quad
      )
      
      n_mod <- length(modelos)
      
      for (i in 1:n_mod) {
        for (j in 1:n_mod) {
          
          m1 <- modelos[[i]]
          m2 <- modelos[[j]]
          
          distancias_df <- distancias_df %>% add_row(
            nbasis = nb,
            modelo1 = names(modelos)[i],
            modelo2 = names(modelos)[j],
            dist_alpha = abs(m1$alpha - m2$alpha),
            dist_beta  = dist_euclidea(m1$beta, m2$beta),
            dist_gamma = dist_euclidea(as.vector(m1$gamma), as.vector(m2$gamma))
          )
          
        }
      }
      
      resultados <- rbind(
        resultados,
        data.frame(
          escenario = basis$nbasis,
          nbasis = basis$nbasis,
          modelo = "fpca_Upsilon0",
          cuadratico = TRUE,
          AUC = auc(roc(y_valid, y_prob_fpca)),
          cond_psi = cond_psi,
          cond_var_fd = cond_var_fd_train
          
          
        ),
        data.frame(
          escenario = basis$nbasis,
          nbasis = basis$nbasis,
          modelo = "fpca_Upsilon0.2_equivalente",
          cuadratico = TRUE,
          AUC = auc(roc(y_valid, y_prob_fpca_equiv)),
          cond_psi = cond_psi,
          cond_var_fd = cond_var_fd_train
          
          
        ),
        data.frame(
          escenario = basis$nbasis,
          nbasis = basis$nbasis,
          modelo = "pca_Upsilon1",
          cuadratico = TRUE,
          AUC = auc(roc(y_valid, y_prob_pca)),
          cond_psi = cond_psi,
          cond_var_fd = cond_var_fd_train
          
        )
      )
    }
  }
  
  return(list(resultados = resultados, distancias = distancias_df))
}


#################################### funciones auxiliares para exportar graficos ############################

#### Exportar gráficos para tesis (sin títulos) ####
  exportar_plot_beta_full <- function(t_grid, beta_true = NULL, beta_est, 
                                      title_suffix = "", output_filename = "beta_plot.png") {
    library(ggplot2)
    library(reshape2)
    
    # Dataframe base
    df_beta <- data.frame(
      t        = t_grid,
      Estimado = beta_est
    )
    
    if (!is.null(beta_true)) {
      df_beta$Original <- beta_true
    }
    
    # Formato largo
    df_beta_long <- reshape2::melt(
      df_beta, id.vars = "t",
      variable.name = "Modelo",
      value.name   = "beta"
    )
    
    # Plot
    p <- ggplot(df_beta_long, aes(x = t, y = beta, color = Modelo, linetype = Modelo)) +
      geom_line(size = 1) +
      scale_color_manual(values = c("Original" = "black", "Estimado" = "darkred")) +
      scale_linetype_manual(values = c("Original" = "solid", "Estimado" = "solid")) +
      labs(x = "t", y = expression(beta(t))) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_blank(),
        
        # ---- NÚMEROS Y TÍTULOS DE EJES EN NEGRITA ----
        axis.text  = element_text(size = 16, colour = "black", face = "bold"),
        axis.title = element_text(size = 16, colour = "black", face = "bold")
      )
    
    ggsave(filename = output_filename, plot = p, width = 6, height = 4, dpi = 300)
    message(paste("Gráfico β(t) guardado como:", output_filename))
  }
  

exportar_gamma_contour_full <- function(t_grid, gamma, title_suffix = "", 
                                        output_filename = "gamma_contour.png") {
  library(ggplot2)
  library(reshape2)
  library(viridis)
  library(grid)  # para unit()
  
  prep_gamma <- function(mat, tipo) {
    df <- reshape2::melt(mat)
    colnames(df) <- c("s_idx", "t_idx", "gamma")
    df$Tipo <- tipo
    df$s    <- t_grid[df$s_idx]
    df$t    <- t_grid[df$t_idx]
    return(df)
  }
  
  if (is.list(gamma)) {
    df_gamma <- do.call(rbind,
                        lapply(names(gamma), function(tipo) prep_gamma(gamma[[tipo]], tipo)))
  } else {
    df_gamma <- prep_gamma(gamma, "gamma")
  }
  
  p <- ggplot(df_gamma, aes(x = t, y = s, z = gamma)) +
    geom_contour_filled(aes(fill = ..level..), bins = 20) +
    facet_wrap(~ Tipo, nrow = 1) +
    scale_fill_viridis_d(option = "plasma", name = expression(gamma(s, t))) +
    labs(x = "t", y = "s") +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 12),
      
      axis.text  = element_text(size = 17, colour = "black", face = "bold"),
      axis.title = element_text(size = 16, colour = "black", face = "bold"),
      
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8),
      
      legend.key.height = unit(0.3, "cm"),
      legend.key.width  = unit(0.3, "cm"),
      legend.spacing.y  = unit(0.2, "cm")
    )
  
  ggsave(filename = output_filename, plot = p, width = 8, height = 4, dpi = 300)
  message(paste("Gráfico γ(s,t) guardado como:", output_filename))
}


exportar_plot_gamma_3d_full <- function(
    t_grid, gamma, title_suffix = "", output_filename = "gamma_3d.png",
    impute_na = FALSE, byrow_guess = TRUE, n_max = 150, palette_name = "Blues") {
  
  library(RColorBrewer)  # para paletas predefinidas
  
  # Parámetros visuales
  theta   <- 35
  phi_ang <- 25
  expand  <- 0.8
  
  # --- Validaciones ---
  if (is.null(t_grid) || length(t_grid) < 2) stop("t_grid debe ser vector numérico >= 2.")
  if (!is.numeric(t_grid)) stop("t_grid debe ser numérico.")
  
  gamma_mat <- if (is.data.frame(gamma)) as.matrix(gamma) else gamma
  if (!is.matrix(gamma_mat)) stop("'gamma' debe ser matriz.")
  if (!is.numeric(gamma_mat)) gamma_mat <- apply(gamma_mat, 2, as.numeric)
  
  n_t <- length(t_grid)
  if (nrow(gamma_mat) != n_t || ncol(gamma_mat) != n_t) {
    stop("Dimensiones de 'gamma' (", nrow(gamma_mat), "x", ncol(gamma_mat),
         ") no coinciden con length(t_grid) = ", n_t)
  }
  
  # --- Submuestreo para visualización ---
  if (n_t > n_max) {
    idx <- seq(1, n_t, length.out = n_max)
    gamma_mat <- gamma_mat[idx, idx]
    t_grid <- t_grid[idx]
  }
  
  # --- Colores tipo RdPu ---
  zmin <- min(gamma_mat, na.rm = TRUE)
  zmax <- max(gamma_mat, na.rm = TRUE)
  ncol <- 100
  colfunc <- colorRampPalette(rev(brewer.pal(9, palette_name)))
  col_z <- colfunc(ncol)[as.numeric(cut(gamma_mat, breaks = ncol))]
  
  # --- Configuración del dispositivo ---
  png(filename = output_filename, width = 1000, height = 850, res = 150, bg = "white")
  old_par <- par(no.readonly = TRUE)
  on.exit({ par(old_par); dev.off() }, add = TRUE)
  
  par(mar = c(1.5, 1.5, 1.2, 0.8), oma = c(0, 0, 0, 0), mgp = c(1.3, 0.3, 0))
  
  # --- Gráfico 3D ---
  persp(
    x = t_grid, y = t_grid, z = gamma_mat,
    theta = theta, phi = phi_ang, expand = expand,
    col = col_z, border = NA,
    box = TRUE, zlim = c(zmin, zmax),
    ticktype = "detailed",
    xlab = "t", ylab = "s", zlab = expression(gamma(s,t)),
    
    # --- Números de los ejes ---
    cex.axis = 1.4,
    col.axis = "black",
    font.axis = 2,     # <<--- NEGRITA
    
    # --- Títulos de ejes ---
    cex.lab = 1.6,
    col.lab = "black",
    font.lab = 2       # <<--- NEGRITA
  )
  
  message(paste("Gráfico 3D guardado como:", output_filename))
}


########################################## funciones auxiliares para operar matrices #########################


# ---- Función auxiliar: vech inverso (reconstruye matriz simétrica) ----
vech_inv <- function(v, p) {
  if (p == 1) {
    if (length(v) != 1)
      stop("Para p = 1, el vector v debe tener longitud 1.")
    return(matrix(v, nrow = 1, ncol = 1))
  }
  
  M <- matrix(0, nrow = p, ncol = p)
  M[lower.tri(M, diag = TRUE)] <- v
  M <- M + t(M) - diag(diag(M))
  return(M)
}

# ---- Función auxiliar: matriz M (1 en diag, 2 en triángulo inferior, 0 arriba) ----
M_matrix <- function(p) {
  M <- matrix(0, nrow = p, ncol = p)
  M[lower.tri(M)] <- 2
  diag(M) <- 1
  return(M)
}

M_matrix_inv <- function(p) {
  M <- matrix(0, nrow = p, ncol = p)
  M[lower.tri(M)] <- 1/2
  diag(M) <- 1
  return(M)
}


#vech^-1(X cdot M^-1)
reconstruir_matriz_simetrica_div2 <- function(vec) {
  # Calcular dimensión
  p <- (sqrt(8 * length(vec) + 1) - 1) / 2
  
  if (p != floor(p)) stop("La longitud del vector no corresponde a un triángulo inferior válido.")
  p <- as.integer(p)
  
  # Crear matriz vacía
  M <- matrix(0, nrow = p, ncol = p)
  
  # Rellenar triángulo inferior
  M[lower.tri(M, diag = TRUE)] <- vec
  
  # Hacerla simétrica
  M[upper.tri(M)] <- t(M)[upper.tri(M)]
  
  # Dividir entre 2 los valores fuera de la diagonal
  M[lower.tri(M)] <- M[lower.tri(M)] / 2
  M[upper.tri(M)] <- M[upper.tri(M)] / 2
  
  return(M)
}


#Q(X) con X matriz
generar_matriz_cuadratica <- function(X) {
  X <- as.matrix(X)
  p <- ncol(X)
  nombres <- colnames(X)
  
  # Lista para acumular resultados
  lista_quad <- list()
  nombres_quad <- c()
  
  for (i in 1:p) {
    for (j in i:p) {
      lista_quad[[length(lista_quad) + 1]] <- X[, i] * X[, j]
      nombres_quad <- c(nombres_quad, paste0(nombres[i], "_x_", nombres[j]))
    }
  }
  
  # Convertir lista a matriz
  X_quad <- do.call(cbind, lista_quad)
  colnames(X_quad) <- nombres_quad
  
  return(X_quad)
}



############################################## GD-penalizado solver general ###################################
gradient_descent_penalized_quad <- function(
    Z_lin, Z_quad, y,
    transformada_beta_fun,
    transformada_gamma_fun,
    H_beta_fun,
    H_gamma_fun,
    modelo_quad = TRUE,
    step_gradient = 0.01,
    iterations = 1000,
    batch_size = 100,
    tol = 1e-6,
    lambda_lin = 0,
    lambda_quad = 0,
    penalty_beta_reduced = diag(ncol(Z_lin)),
    penalty_gamma_reduced = diag(ncol(Z_quad)),
    beta = c(),
    gamma = matrix(0, nrow = 1, ncol = 1),
    seed = 123,
    verbose = TRUE,
    V_lin = NULL,
    V_quad = NULL,
    psi = NULL
) {
  set.seed(seed)
  n <- nrow(Z_lin)
  K_lin <- ncol(Z_lin)
  K_quad <- ncol(Z_quad)
  
  alpha <- 1e-8
  
  # Inicialización flexible según si V_lin / V_quad son NULL
  if (is.null(beta) || length(beta) != K_lin) {
    UpsBeta <- rep(1e-8, K_lin)
  } else {
    args_beta <- list(beta)
    if (!is.null(V_lin)) args_beta <- c(args_beta, list(V_lin))
    if (!is.null(psi)) args_beta <- c(args_beta, list(psi))
    UpsBeta <- do.call(transformada_beta_fun, args_beta)
  }
  
  if (is.null(gamma) || !all(dim(gamma) == c(K_lin, K_lin))) {
    UpsGamma <- rep(1e-8, K_quad)
  } else {
    args_gamma <- list(gamma)
    if (!is.null(V_quad))        args_gamma <- c(args_gamma, list(V_quad))
    if (!is.null(psi))  args_gamma <- c(args_gamma, list(psi))
    UpsGamma <- do.call(transformada_gamma_fun, args_gamma)
  }
  
  cost_history <- numeric(iterations)
  
  for (i in 1:iterations) {
    idx <- sample(1:n, batch_size)
    Zb_lin <- Z_lin[idx, , drop = FALSE]
    Zb_quad <- Z_quad[idx, , drop = FALSE]
    yb <- y[idx]
    nb <- length(idx)
    
    lin_pred <- Zb_lin %*% UpsBeta
    quad_pred <- if (modelo_quad) Zb_quad %*% UpsGamma else 0
    x <- pmax(pmin(alpha + lin_pred + quad_pred, 20), -20)
    y_pred <- 1 / (1 + exp(-x))
    residuos <- yb - y_pred
    
    # --- Gradiente alpha ---
    grad_alpha <- -(1/nb) * sum(residuos)
    
    # --- Gradiente beta ---
    args_Hbeta <- list(UpsBeta)
    if (!is.null(V_lin)) args_Hbeta <- c(args_Hbeta, list(V_lin))
    if (!is.null(psi))  args_Hbeta <- c(args_Hbeta, list(psi))
    Hbeta_val <- do.call(H_beta_fun, args_Hbeta)
    
    penal_lin <- lambda_lin * 2 * penalty_beta_reduced %*% Hbeta_val
    args_transform_beta <- list(penal_lin)
    if (!is.null(V_lin)) args_transform_beta <- c(args_transform_beta, list(V_lin))
    if (!is.null(psi)) args_transform_beta <- c(args_transform_beta, list(psi))
    grad_beta <- -(1/nb) * t(Zb_lin) %*% residuos + do.call(transformada_beta_fun, args_transform_beta)
    
    # --- Gradiente gamma ---
    if (modelo_quad) {
      args_Hgamma <- list(UpsGamma)
      if (!is.null(V_quad)) args_Hgamma <- c(args_Hgamma, list(V_quad))
      if (!is.null(psi)) args_Hgamma <- c(args_Hgamma, list(psi) )
      Hgamma_val <- do.call(H_gamma_fun, args_Hgamma)
      
      penal_quad <- lambda_quad * 2 * penalty_gamma_reduced %*% Hgamma_val %*% t(penalty_gamma_reduced)
      args_transform_gamma <- list(penal_quad)
      if (!is.null(V_quad)) args_transform_gamma <- c(args_transform_gamma, list(V_quad))
      if (!is.null(psi)) args_transform_gamma <- c(args_transform_gamma, list(psi))
      grad_gamma <- -(1/nb) * t(Zb_quad) %*% residuos + do.call(transformada_gamma_fun, args_transform_gamma)
    }
    
    # --- Actualización ---
    alpha <- alpha - step_gradient * grad_alpha
    UpsBeta <- UpsBeta - step_gradient * grad_beta
    if (modelo_quad) UpsGamma <- UpsGamma - step_gradient * grad_gamma
    
    # --- Evaluación costo ---
    x_full <- pmax(pmin(alpha + Z_lin %*% UpsBeta + if (modelo_quad) Z_quad %*% UpsGamma else 0, 20), -20)
    y_pred_full <- 1 / (1 + exp(-x_full))
    eps <- 1e-8
    y_pred_full <- pmin(pmax(y_pred_full, eps), 1 - eps)
    cost <- -mean(y * log(y_pred_full) + (1 - y) * log(1 - y_pred_full))
    cost_history[i] <- cost
    
    if (i > 1 && abs(cost_history[i] - cost_history[i - 1]) < tol) {
      if (verbose) cat("Convergencia alcanzada en iteración", i, "\n")
      break
    }
    
    if (verbose && (i %% 10 == 0)) {
      cat("Iter:", i, "Costo:", round(cost, 6), "\n")
    }
  }
  
  list(
    alpha = alpha,
    UpsBeta = UpsBeta,
    UpsGamma = UpsGamma,
    cost_history = cost_history
  )
}


################################################# naive ########################################


H_beta_naive <- function(beta, return_matrix = FALSE){
  return(beta)
}
transformada_beta_naive <- function(beta){
  return(beta)
}

H_gamma_naive <- function(Upsilon){
  Upsilon <- as.vector(Upsilon)
  p <- (-1 + sqrt(1 + 8 * length(Upsilon))) / 2
  # construir (M/2) 
  M <- M_matrix_inv(p)
  M_vec <- M[lower.tri(M, diag = TRUE)] 
  
  # reconstruir matriz simétrica
  Gamma <- vech_inv(Upsilon * M_vec, p)
  
  # aplicar la transformación
  return(Gamma)
}

transformada_gamma_naive <- function(gamma, return_matrix = FALSE){
  # construir M: 1 en diagonal, 2 en triángulo inferior, 0 en triángulo superior
  p <- nrow(gamma)
  M <- matrix(0, nrow = p, ncol = p)
  diag(M) <- 1
  M[lower.tri(M)] <- 2
  
  # multiplicación elemento a elemento
  B <- gamma * M
  
  if (return_matrix) {
    return(B)   # si se quiere la matriz completa
  }
  
  # vech: tomar triángulo inferior por columnas (incluye diagonal)
  v <- B[lower.tri(B, diag = TRUE)]
  return(as.numeric(v))
}


Upsilon_naive <- function(             fd_centered,
                                       y,
                                       beta,
                                       gamma,
                                       alpha,
                                       step_gradient,
                                       iterations,
                                       basis,
                                       LdPenalization,
                                       lambda_lin = 1e-12, 
                                       lambda_quad = 1e-12, 
                                       modelo_quad = TRUE,
                                       tol = 1e-6,
                                       batch_size = 100,
                                       verbose = TRUE) {   
  
  
  
  cost_history <- numeric(iterations)
  psi <- inprod(basis, basis)
  
  # 1. Matriz A y A_psi
  A <- t(fd_centered$coefs)
  A_psi <- A %*% psi
  
  
  # 2. Rugosidad
  LdBasis <- getbasispenalty(basis, LdPenalization)
  penalty_beta_reduced <- LdBasis      # K x K
  penalty_beta_reduced <- penalty_beta_reduced / norm(penalty_beta_reduced, type = "F")
  penalty_gamma_reduced <- LdBasis    # K x K
  penalty_gamma_reduced <- penalty_gamma_reduced / norm(penalty_gamma_reduced, type = "F")
  
  res_gradiente <- gradient_descent_penalized_quad(     Z_lin = A_psi,
                                                         Z_quad = generar_matriz_cuadratica(A_psi),
                                                         y = y,
                                                         transformada_beta_fun = transformada_beta_naive,
                                                         transformada_gamma_fun = transformada_gamma_naive,
                                                         H_beta_fun = H_beta_naive,
                                                         H_gamma_fun = H_gamma_naive,
                                                         modelo_quad = modelo_quad,
                                                         step_gradient = step_gradient,
                                                         tol = tol,
                                                         batch_size = batch_size,
                                                         lambda_lin = lambda_lin,
                                                         lambda_quad = lambda_quad,
                                                         iterations = iterations,
                                                         penalty_beta_reduced = penalty_beta_reduced,
                                                         penalty_gamma_reduced = penalty_gamma_reduced,
                                                         beta = beta,
                                                         gamma = gamma,
                                                        verbose = verbose)
  
  gamma_original  <- H_gamma_naive(res_gradiente$UpsGamma)
  beta_original   <- res_gradiente$UpsBeta 
  alpha_original  <- res_gradiente$alpha
  
  return(list(
    beta_reduced = res_gradiente$UpsBeta,
    gamma_reduced = res_gradiente$UpsGamma,
    beta = as.vector(beta_original),
    gamma = gamma_original,
    alpha = alpha_original,
    cost_history = res_gradiente$cost_history
  ))
}
################################################ SOBRE BASE FPCA UPSILON 0 ###########################


H_beta_0 <- function(Upsilon) {
    return(Upsilon)
}

# ---- Operador H_gamma^{(1)}(Upsilon) = V^{(1)} * vech^{-1}(Upsilon * M/2) * V^{(1)T} ----
H_gamma_0 <- function(Upsilon) {
  Upsilon <- as.vector(Upsilon)
  p <- (-1 + sqrt(1 + 8 * length(Upsilon))) / 2
  # construir (M/2) 
  M <- M_matrix_inv(p)
  M_vec <- M[lower.tri(M, diag = TRUE)] 
  
  # reconstruir matriz simétrica
  Gamma <- vech_inv(Upsilon * M_vec, p)
  
  # aplicar la transformación
  return(Gamma)
}

transformada_beta_fpca <- function(beta) {
  beta <- as.matrix(beta)
  return(beta)
}

# vech_transform: devuelve vech(gamma * M)
transformada_gamma_fpca <- function(gamma, return_matrix = FALSE) {
  # construir M: 1 en diagonal, 2 en triángulo inferior, 0 en triángulo superior
  p <- ncol(gamma)
  M <- matrix(0, nrow = p, ncol = p)
  diag(M) <- 1
  M[lower.tri(M)] <- 2
  
  # multiplicación elemento a elemento
  B <- gamma * M
  
  if (return_matrix) {
    return(B)   # si se quiere la matriz completa
  }
  
  # vech: tomar triángulo inferior por columnas (incluye diagonal)
  v <- B[lower.tri(B, diag = TRUE)]
  return(as.numeric(v))
}



fpca_Upsilon0 <- function( fd_centered,
                           y,
                           beta,
                           gamma,
                           alpha,
                           step_gradient,
                           iterations,
                           K = NULL,                    # <<< NEW ARGUMENT
                           var_threshold = 0.95,
                           basis,
                           LdPenalization,
                           lambda_lin = 1e-12, 
                           lambda_quad = 1e-12, 
                           modelo_quad = TRUE,
                           tol = 1e-6,
                           batch_size = 100,
                           scale = TRUE,
                           verbose = TRUE)
{
  
  cost_history <- numeric(iterations)
  
  n <- nrow(fd_centered$coef)
  
  # --- Primera fpca para obtener proporciones ---
  fpca_temp <- pca.fd(fd_centered, nharm = fd_centered$basis$nbasis, centerfns = TRUE)
  cum_var <- cumsum(fpca_temp$varprop)
  
  # --- elegir K si no viene ---
  if (is.null(K)) {
    K <- which(cum_var >= var_threshold)[1]
    
    if (is.na(K)) {
      K <- length(cum_var)
      warning(sprintf(
        "No se alcanzó el umbral de varianza %.2f, usando todas las %d componentes.",
        var_threshold, K
      ))
    }
  }
  
  # sanity check
  if (K <= 0) stop("K debe ser >= 1")
  if (K > ncol(fpca_temp$scores)) stop("K mayor que número de componentes válidas FPCA")
  
  # --- fpca truncada final ---
  fpca_result <- pca.fd(fd_centered, nharm = K, centerfns = TRUE)
  
  Xi_hat <- fpca_result$scores
  phi_hat <- fpca_result$harmonics
  
  Z_lin <- Xi_hat[, 1:K, drop=FALSE]
  Z_quad <- generar_matriz_cuadratica(Z_lin)
  
  # Rugosidad
  basis_phi <- phi_hat$basis
  LdBasis <- getbasispenalty(basis_phi, LdPenalization)
  LdBasisNorm <- LdBasis / norm(LdBasis, "F")
  R <- t(phi_hat$coef) %*% LdBasisNorm %*% phi_hat$coef
  
  res_gradiente <- gradient_descent_penalized_quad(
    Z_lin = Z_lin,
    Z_quad = Z_quad,
    y = y,
    transformada_beta_fun = transformada_beta_fpca,
    transformada_gamma_fun = transformada_gamma_fpca,
    H_beta_fun = H_beta_0,
    H_gamma_fun = H_gamma_0,
    modelo_quad = modelo_quad,
    step_gradient = step_gradient,
    tol = tol,
    batch_size = batch_size,
    lambda_lin = lambda_lin,
    lambda_quad = lambda_quad,
    iterations = iterations,
    penalty_beta_reduced = R,
    penalty_gamma_reduced = R,
    beta = beta,
    gamma = gamma,
    verbose = verbose)
  
  gamma_original  <- H_gamma_0(res_gradiente$UpsGamma)
  beta_original   <- H_beta_0(res_gradiente$UpsBeta) 
  alpha_original  <- res_gradiente$alpha
  
  return(list(
    beta_reduced = res_gradiente$UpsBeta,
    gamma_reduced = res_gradiente$UpsGamma,
    beta = as.vector(beta_original),
    gamma = gamma_original,
    alpha = alpha_original,
    cost_history = res_gradiente$cost_history,
    K = K,
    fpca_model = fpca_result
  ))
}

################################################ FPCA Equivalencia UPSILON 0 .2 ###########################
Psi_inv_half <- function(Psi) {
  eig <- eigen(Psi, symmetric = TRUE)
  V <- eig$vectors
  L <- eig$values
  return(V %*% diag(1 / sqrt(L)) %*% t(V))

}
Psi_half <- function(Psi) {
  eig <- eigen(Psi, symmetric = TRUE)
  V <- eig$vectors
  L <- eig$values
  return(V %*% diag(sqrt(L)) %*% t(V))
}


# ---- Operador H_beta^{(0.2)}(Upsilon) = psi^{-1/2}V^{(0.2)} * Upsilon ----
H_beta_0_equiv <- function(Upsilon,V,psi) {
  V <- as.matrix(V)
  Upsilon <- as.vector(Upsilon)
  psi_half_inv <- Psi_inv_half(psi)
  return(psi_half_inv %*% V %*% Upsilon)
}

# ---- Operador H_gamma^{(0.2)}(Upsilon) = psi^{-1/2}V^{(0.2)} * vech^{-1}(Upsilon * M/2) * V^{(0.2)T}psi^{-1/2} ----
H_gamma_0_equiv  <- function(Upsilon, V, psi) {
  V <- as.matrix(V)
  Upsilon <- as.vector(Upsilon)
  
  # tamaño implícito de matriz simétrica
  p <- ncol(V)
  len_expected <- p * (p + 1) / 2
  
  if (length(Upsilon) != len_expected) {
    stop(sprintf("Longitud de Upsilon incorrecta: esperada %d (para p=%d).", len_expected, p))
  }
  
  # construir (M/2) 
  M <- M_matrix_inv(p)
  M_vec <- M[lower.tri(M, diag = TRUE)] 
  
  # reconstruir matriz simétrica
  Gamma <- vech_inv(Upsilon * M_vec, p)
  
  psi_half_inv <- Psi_inv_half(psi)
  
  # aplicar la transformación
  return(psi_half_inv%*%V %*% Gamma %*% t(V)%*%psi_half_inv)
}


# vT_beta: calcula psi^{1/2}%*% t(V) %*% beta
transformada_beta_fpca_equiv  <- function(beta, V, psi) {
  V <- as.matrix(V)
  beta <- as.matrix(beta)
  
  # chequeo de compatibilidad de dimensiones
  if (nrow(V) != nrow(beta)) {
    stop(sprintf("Dimensiones incompatibles: nrow(V) = %d, nrow(beta) = %d", nrow(V), nrow(beta)))
  }
  
  # producto V^T * beta
  result <- t(V) %*% Psi_half(psi) %*% beta
  return(result)
}



# vech_transform: devuelve vech( (t(V) %*% gamma %*% V)  *  M )
transformada_gamma_fpca_equiv <- function(gamma, V, psi, return_matrix = FALSE) {
  # coerciones y checks básicos
  V <- as.matrix(V)
  gamma <- as.matrix(gamma)
  if (!is.numeric(V) || !is.numeric(gamma)) {
    stop("V y gamma deben ser matrices numéricas.")
  }
  # dimensiones compatibles: gamma debe ser p x p, V debe ser p x p' o p x k tal que t(V) %*% gamma %*% V tenga sentido.
  if (nrow(V) != nrow(gamma)) {
    stop(sprintf("Dimensiones incompatibles: ncol(V) = %d, nrow(gamma) = %d", ncol(V), nrow(gamma)))
  }
  
  psi_half <- Psi_half(psi)
  # calcular A = t(V) %*% gamma %*% V
  A <- t(V) %*% psi_half %*% gamma %*% psi_half %*% V       
  
  # construir M: 1 en diagonal, 2 en triángulo inferior, 0 en triángulo superior
  p <- nrow(A)
  M <- matrix(0, nrow = p, ncol = p)
  diag(M) <- 1
  M[lower.tri(M)] <- 2
  
  # multiplicación elemento a elemento
  B <- A * M
  
  if (return_matrix) {
    return(B)   # si se quiere la matriz completa
  }
  
  # vech: tomar triángulo inferior por columnas (incluye diagonal)
  v <- B[lower.tri(B, diag = TRUE)]
  return(as.numeric(v))
}



fpca_upsilon0_equiv<- function( fd_centered,
                                y,
                                beta,
                                gamma,
                                alpha,
                                step_gradient,
                                iterations,
                                K = NULL,                 # <<--- nuevo parámetro
                                var_threshold = 0.95,
                                basis,
                                LdPenalization,
                                lambda_lin = 1e-12, 
                                lambda_quad = 1e-12, 
                                modelo_quad = TRUE,
                                tol = 1e-6,
                                batch_size = 100,
                                scale = TRUE,
                                verbose = TRUE) 
{   
  cost_history <- numeric(iterations)
  psi <- inprod(basis, basis)
  
  # 1) Matriz A y A_psi_half
  A <- t(fd_centered$coefs)
  A_psi_half <- A %*% Psi_half(psi)
  
  # 2) PCA sobre A_psi_half
  pca_result <- prcomp(A_psi_half, center = TRUE, scale. = scale)
  
  # 3) Determinar K correcto
  if (is.null(K)) {
    var_exp_acum <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
    K <- which(var_exp_acum >= var_threshold)[1]
  }
  
  # sanity check
  if (K > ncol(pca_result$rotation)) {
    stop("Valor de K mayor al número de componentes disponibles.")
  }
  
  # extraer coordenadas
  Z_lin <- pca_result$x[, 1:K, drop=FALSE]
  Z_quad <- generar_matriz_cuadratica(Z_lin)
  
  # 5. Rugosidad
  LdBasis <- getbasispenalty(basis, LdPenalization)
  
  V <- pca_result$rotation[, 1:K, drop=FALSE]
  
  penalty_beta_reduced <- LdBasis / norm(LdBasis, type="F")
  penalty_gamma_reduced <- LdBasis / norm(LdBasis, type="F")
  
  res_gradiente <- gradient_descent_penalized_quad(    
    Z_lin = Z_lin,
    Z_quad = Z_quad,
    y = y,
    transformada_beta_fun = transformada_beta_fpca_equiv,
    transformada_gamma_fun = transformada_gamma_fpca_equiv,
    H_beta_fun = H_beta_0_equiv,
    H_gamma_fun = H_gamma_0_equiv,
    modelo_quad = modelo_quad,
    step_gradient = step_gradient,
    tol = tol,
    batch_size = batch_size,
    lambda_lin = lambda_lin,
    lambda_quad = lambda_quad,
    iterations = iterations,
    penalty_beta_reduced = penalty_beta_reduced,
    penalty_gamma_reduced = penalty_gamma_reduced,
    beta = beta,
    gamma = gamma,
    V_lin = V,
    V_quad = V,
    psi = psi,
    verbose = verbose)
  
  # Ajuste V original desescalado
  mu <- pca_result$center
  s  <- if (scale) pca_result$scale else rep(1, length(mu))
  s_safe <- ifelse(s==0,1,s)
  
  V <- (pca_result$rotation * 1/s_safe)[, 1:K, drop=FALSE]
  
  gamma_original  <- H_gamma_0_equiv(res_gradiente$UpsGamma, V, psi)
  beta_original   <- H_beta_0_equiv(res_gradiente$UpsBeta, V, psi)
  alpha_original  <- res_gradiente$alpha
  
  return(list(
    beta_reduced = res_gradiente$UpsBeta,
    gamma_reduced = res_gradiente$UpsGamma,
    beta = as.vector(beta_original),
    gamma = gamma_original,
    alpha = alpha_original,
    cost_history = res_gradiente$cost_history,
    K = K,
    pca_model = pca_result,
    coords_f_en_base = Psi_half(psi) %*% V
  ))
}


# ---- Operador H_beta^{(0.2)}(Upsilon) = psi^{-1/2}V^{(0.2)} * Upsilon ----
H_beta_0_equiv_normalizado <- function(Upsilon,V, M) {
  V <- as.matrix(V)
  Upsilon <- as.vector(Upsilon)
  Minv <- solve(M)
  return(Minv %*% V %*% Upsilon)
}

# ---- Operador H_gamma^{(0.2)}(Upsilon) = psi^{-1/2}V^{(0.2)} * vech^{-1}(Upsilon * M/2) * V^{(0.2)T}psi^{-1/2} ----
H_gamma_0_equiv_normalizado  <- function(Upsilon, V, M) {
  V <- as.matrix(V)
  Upsilon <- as.vector(Upsilon)
  Minv <- solve(M)
  
  # tamaño implícito de matriz simétrica
  p <- ncol(V)
  len_expected <- p * (p + 1) / 2
  
  if (length(Upsilon) != len_expected) {
    stop(sprintf("Longitud de Upsilon incorrecta: esperada %d (para p=%d).", len_expected, p))
  }
  
  # construir (M/2) 
  M <- M_matrix_inv(p)
  M_vec <- M[lower.tri(M, diag = TRUE)] 
  
  # reconstruir matriz simétrica
  Gamma <- vech_inv(Upsilon * M_vec, p)
  
  # aplicar la transformación
  return(Minv%*%V %*% Gamma %*% t(V)%*%t(Minv))
}


# vT_beta: calcula psi^{1/2}%*% t(V) %*% beta
transformada_beta_fpca_equiv_normalizado  <- function(beta, V, M) {
  V <- as.matrix(V)
  beta <- as.matrix(beta)
  
  # chequeo de compatibilidad de dimensiones
  if (nrow(V) != nrow(beta)) {
    stop(sprintf("Dimensiones incompatibles: nrow(V) = %d, nrow(beta) = %d", nrow(V), nrow(beta)))
  }
  
  # producto V^T * beta
  result <- t(V) %*% M %*% beta
  return(result)
}



# vech_transform: devuelve vech( (t(V) %*% gamma %*% V)  *  M )
transformada_gamma_fpca_equiv_normalizado <- function(gamma, V, M, return_matrix = FALSE) {
  # coerciones y checks básicos
  V <- as.matrix(V)
  gamma <- as.matrix(gamma)
  if (!is.numeric(V) || !is.numeric(gamma)) {
    stop("V y gamma deben ser matrices numéricas.")
  }
  # dimensiones compatibles: gamma debe ser p x p, V debe ser p x p' o p x k tal que t(V) %*% gamma %*% V tenga sentido.
  if (nrow(V) != nrow(gamma)) {
    stop(sprintf("Dimensiones incompatibles: ncol(V) = %d, nrow(gamma) = %d", ncol(V), nrow(gamma)))
  }
  
  A <- t(V) %*% M %*% gamma %*% t(M) %*% V       
  
  # construir M: 1 en diagonal, 2 en triángulo inferior, 0 en triángulo superior
  p <- nrow(A)
  M <- matrix(0, nrow = p, ncol = p)
  diag(M) <- 1
  M[lower.tri(M)] <- 2
  
  # multiplicación elemento a elemento
  B <- A * M
  
  if (return_matrix) {
    return(B)   # si se quiere la matriz completa
  }
  
  # vech: tomar triángulo inferior por columnas (incluye diagonal)
  v <- B[lower.tri(B, diag = TRUE)]
  return(as.numeric(v))
}



#gradient_descent_penalized_PCA -> pca_coef_Upsilon1
fpca_upsilon0_equiv_normalizado <- function( fd_centered,
                                             y,
                                             beta,
                                             gamma,
                                             alpha,
                                             step_gradient,
                                             iterations,
                                             K = NULL,                 # <<< NEW
                                             var_threshold = 0.95,
                                             basis,
                                             LdPenalization,
                                             lambda_lin = 1e-12, 
                                             lambda_quad = 1e-12, 
                                             modelo_quad = TRUE,
                                             tol = 1e-6,
                                             batch_size = 100,
                                             scale = TRUE,
                                             verbose = TRUE) 
{   
  
  cost_history <- numeric(iterations)
  psi <- inprod(basis, basis)
  
  M <- chol(psi)  
  Minv <- solve(M)
  
  psi <- Minv %*% psi %*% t(Minv)
  
  A <- t(fd_centered$coefs)
  A_psi_half <- A %*% t(M) %*% Psi_half(psi)
  
  pca_result <- prcomp(A_psi_half, center = TRUE, scale. = scale)
  
  ## seleccionar K si no viene
  if (is.null(K)) {
    var_exp_acum <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
    K <- which(var_exp_acum >= var_threshold)[1]
    
    if (is.na(K)) {
      K <- ncol(pca_result$rotation)
      warning(sprintf("No se alcanzó %.2f de varianza — usando K = %d", 
                      var_threshold, K))
    }
  }
  
  if (K <= 0) stop("K debe ser >=1")
  if (K > ncol(pca_result$rotation)) stop("K mayor al número posible de componentes")
  
  Z_lin <- pca_result$x[, 1:K, drop=FALSE]
  Z_quad <- generar_matriz_cuadratica(Z_lin)
  
  LdBasis <- getbasispenalty(basis, LdPenalization)
  
  V <- pca_result$rotation[, 1:K, drop=FALSE]
  
  penalty_beta_reduced <- LdBasis / norm(LdBasis,"F")
  penalty_gamma_reduced <- LdBasis / norm(LdBasis,"F")
  
  res_gradiente <- gradient_descent_penalized_quad(
    Z_lin = Z_lin,
    Z_quad = Z_quad,
    y = y,
    transformada_beta_fun = transformada_beta_fpca_equiv_normalizado,
    transformada_gamma_fun = transformada_gamma_fpca_equiv_normalizado,
    H_beta_fun = H_beta_0_equiv_normalizado,
    H_gamma_fun = H_gamma_0_equiv_normalizado,
    modelo_quad = modelo_quad,
    step_gradient = step_gradient,
    tol = tol,
    batch_size = batch_size,
    lambda_lin = lambda_lin,
    lambda_quad = lambda_quad,
    iterations = iterations,
    penalty_beta_reduced = penalty_beta_reduced,
    penalty_gamma_reduced = penalty_gamma_reduced,
    beta = beta,
    gamma = gamma,
    V_lin = V,
    V_quad = V,
    psi = M,
    verbose = verbose)
  
  mu <- pca_result$center
  s  <- if (scale) pca_result$scale else rep(1,length(mu))
  s_safe <- ifelse(s==0,1,s)
  
  V <- (pca_result$rotation * 1/s_safe)[,1:K,drop=FALSE]
  
  gamma_original <- H_gamma_0_equiv(res_gradiente$UpsGamma, V, M)
  beta_original  <- H_beta_0_equiv(res_gradiente$UpsBeta, V, M)
  
  alpha_original <- res_gradiente$alpha
  
  return(list(
    beta_reduced = res_gradiente$UpsBeta,
    gamma_reduced = res_gradiente$UpsGamma,
    beta = as.vector(beta_original),
    gamma = gamma_original,
    alpha = alpha_original,
    cost_history = res_gradiente$cost_history,
    K = K,
    pca_model = pca_result,
    coords_f_en_base = Psi_half(psi) %*% V
  ))
}


################################################ PCA COEF UPSILON 1 ########################################

# ---- Operador H_beta^{(1)}(Upsilon) = V^{(1)} * Upsilon ----
H_beta_1 <- function(Upsilon,V) {
  V <- as.matrix(V)
  Upsilon <- as.vector(Upsilon)
  
  return(V %*% Upsilon)
}

# ---- Operador H_gamma^{(1)}(Upsilon) = V^{(1)} * vech^{-1}(Upsilon * M/2) * V^{(1)T} ----
H_gamma_1 <- function(Upsilon, V) {
  V <- as.matrix(V)
  Upsilon <- as.vector(Upsilon)
  
  # tamaño implícito de matriz simétrica
  p <- ncol(V)
  len_expected <- p * (p + 1) / 2
  
  if (length(Upsilon) != len_expected) {
    stop(sprintf("Longitud de Upsilon incorrecta: esperada %d (para p=%d).", len_expected, p))
  }
  
  # construir (M/2) 
  M <- M_matrix_inv(p)
  M_vec <- M[lower.tri(M, diag = TRUE)] 
  
  # reconstruir matriz simétrica
  Gamma <- vech_inv(Upsilon * M_vec, p)
  
  # aplicar la transformación
  return(V %*% Gamma %*% t(V))
}


# vT_beta: calcula t(V) %*% beta
transformada_beta_pca <- function(beta, V) {
  V <- as.matrix(V)
  beta <- as.matrix(beta)
  
  # chequeo de compatibilidad de dimensiones
  if (nrow(V) != nrow(beta)) {
    stop(sprintf("Dimensiones incompatibles: nrow(V) = %d, nrow(beta) = %d", nrow(V), nrow(beta)))
  }
  
  # producto V^T * beta
  result <- t(V) %*% beta
  return(result)
}



# vech_transform: devuelve vech( (t(V) %*% gamma %*% V)  *  M )
transformada_gamma_pca <- function(gamma, V, return_matrix = FALSE) {
  # coerciones y checks básicos
  V <- as.matrix(V)
  gamma <- as.matrix(gamma)
  if (!is.numeric(V) || !is.numeric(gamma)) {
    stop("V y gamma deben ser matrices numéricas.")
  }
  # dimensiones compatibles: gamma debe ser p x p, V debe ser p x p' o p x k tal que t(V) %*% gamma %*% V tenga sentido.
  if (nrow(V) != nrow(gamma)) {
    stop(sprintf("Dimensiones incompatibles: ncol(V) = %d, nrow(gamma) = %d", ncol(V), nrow(gamma)))
  }
  
  # calcular A = t(V) %*% gamma %*% V
  A <- t(V) %*% gamma %*% V       # resultado: p x p if V is p x k? -> see check above
  
  # construir M: 1 en diagonal, 2 en triángulo inferior, 0 en triángulo superior
  p <- nrow(A)
  M <- matrix(0, nrow = p, ncol = p)
  diag(M) <- 1
  M[lower.tri(M)] <- 2
  
  # multiplicación elemento a elemento
  B <- A * M
  
  if (return_matrix) {
    return(B)   # si se quiere la matriz completa
  }
  
  # vech: tomar triángulo inferior por columnas (incluye diagonal)
  v <- B[lower.tri(B, diag = TRUE)]
  return(as.numeric(v))
}



pca_coef_Upsilon1 <- function( fd_centered,
                               y,
                               beta,
                               gamma,
                               alpha,
                               step_gradient,
                               iterations,
                               K = NULL,                 # <<< NEW
                               var_threshold = 0.95,
                               basis,
                               LdPenalization,
                               lambda_lin = 1e-12, 
                               lambda_quad = 1e-12, 
                               modelo_quad = TRUE,
                               tol = 1e-6,
                               batch_size = 100,
                               scale = TRUE,
                               verbose = TRUE) 
{   
  
  cost_history <- numeric(iterations)
  psi <- inprod(basis, basis)
  
  # 1. A y A_psi
  A <- t(fd_centered$coefs)
  A_psi <- A %*% psi
  
  # 2. PCA multivariada
  pca_result <- prcomp(A_psi, center = TRUE, scale. = scale)
  
  # 3. Selección de K
  if (is.null(K)) {
    var_exp_acum <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
    K <- which(var_exp_acum >= var_threshold)[1]
    
    if (is.na(K)) {
      K <- ncol(pca_result$rotation)
      warning(sprintf(
        "No se alcanzó el %.2f de varianza — usando K = %d componentes",
        var_threshold, K
      ))
    }
  }
  
  if (K <= 0) stop("K debe ser >= 1.")
  if (K > ncol(pca_result$rotation)) stop("K mayor a la cantidad de componentes disponibles.")
  
  Z_lin <- pca_result$x[, 1:K, drop=FALSE]
  Z_quad <- generar_matriz_cuadratica(Z_lin)
  
  # 5. Rugosidad
  LdBasis <- getbasispenalty(basis, LdPenalization)
  
  V <- pca_result$rotation[, 1:K, drop=FALSE]
  
  penalty_beta_reduced <- LdBasis / norm(LdBasis, "F")
  penalty_gamma_reduced <- LdBasis / norm(LdBasis, "F")
  
  res_gradiente <- gradient_descent_penalized_quad(
    Z_lin = Z_lin,
    Z_quad = Z_quad,
    y = y,
    transformada_beta_fun = transformada_beta_pca,
    transformada_gamma_fun = transformada_gamma_pca,
    H_beta_fun = H_beta_1,
    H_gamma_fun = H_gamma_1,
    modelo_quad = modelo_quad,
    step_gradient = step_gradient,
    tol = tol,
    batch_size = batch_size,
    lambda_lin = lambda_lin,
    lambda_quad = lambda_quad,
    iterations = iterations,
    penalty_beta_reduced = penalty_beta_reduced,
    penalty_gamma_reduced = penalty_gamma_reduced,
    beta = beta,
    gamma = gamma,
    V_lin = V,
    V_quad = V,
    verbose = verbose)
  
  # --- Deshacer scaling ---
  
  mu <- pca_result$center
  s  <- if (scale) pca_result$scale else rep(1, length(mu))
  
  s_safe <- ifelse(s==0,1,s)
  
  V <- (pca_result$rotation * 1/s_safe)[, 1:K, drop=FALSE]
  
  gamma_original <- H_gamma_1(res_gradiente$UpsGamma, V)
  beta_original  <- H_beta_1(res_gradiente$UpsBeta, V)
  alpha_original <- res_gradiente$alpha
  
  return(list(
    beta_reduced = res_gradiente$UpsBeta,
    gamma_reduced = res_gradiente$UpsGamma,
    beta = as.vector(beta_original),
    gamma = gamma_original,
    alpha = alpha_original,
    cost_history = res_gradiente$cost_history,
    K = K,
    pca_model = pca_result
  ))
}

# ---- Operador H_beta^{(1)}(Upsilon) = V^{(1)} * Upsilon ----
H_beta_1_normalizado <- function(Upsilon,V, M) {
  V <- as.matrix(V)
  Upsilon <- as.vector(Upsilon)
  Minv <- solve(M)
  
  return(Minv %*% V %*% Upsilon)
}

# ---- Operador H_gamma^{(1)}(Upsilon) = V^{(1)} * vech^{-1}(Upsilon * M/2) * V^{(1)T} ----
H_gamma_1_normalizado <- function(Upsilon, V, M) {
  V <- as.matrix(V)
  Upsilon <- as.vector(Upsilon)
  Minv <- solve(M)
  # tamaño implícito de matriz simétrica
  p <- ncol(V)
  len_expected <- p * (p + 1) / 2
  
  if (length(Upsilon) != len_expected) {
    stop(sprintf("Longitud de Upsilon incorrecta: esperada %d (para p=%d).", len_expected, p))
  }
  
  # construir (M/2) 
  M <- M_matrix_inv(p)
  M_vec <- M[lower.tri(M, diag = TRUE)] 
  
  # reconstruir matriz simétrica
  Gamma <- vech_inv(Upsilon * M_vec, p)
  
  # aplicar la transformación
  return( Minv %*% V %*% Gamma %*% t(V) %*% t(Minv))
}


# vT_beta: calcula t(V) %*% beta
transformada_beta_pca_normalizado <- function(beta, V, M) {
  V <- as.matrix(V)
  beta <- as.matrix(beta)
  
  # chequeo de compatibilidad de dimensiones
  if (nrow(V) != nrow(beta)) {
    stop(sprintf("Dimensiones incompatibles: nrow(V) = %d, nrow(beta) = %d", nrow(V), nrow(beta)))
  }
  
  # producto V^T * beta
  result <- t(V) %*% M %*% beta
  return(result)
}



# vech_transform: devuelve vech( (t(V) %*% gamma %*% V)  *  M )
transformada_gamma_pca_normalizado <- function(gamma, V, M, return_matrix = FALSE) {
  # coerciones y checks básicos
  V <- as.matrix(V)
  gamma <- as.matrix(gamma)
  if (!is.numeric(V) || !is.numeric(gamma)) {
    stop("V y gamma deben ser matrices numéricas.")
  }
  # dimensiones compatibles: gamma debe ser p x p, V debe ser p x p' o p x k tal que t(V) %*% gamma %*% V tenga sentido.
  if (nrow(V) != nrow(gamma)) {
    stop(sprintf("Dimensiones incompatibles: ncol(V) = %d, nrow(gamma) = %d", ncol(V), nrow(gamma)))
  }
  
  # calcular A = t(V) %*% gamma %*% V
  A <- t(V) %*% M %*% gamma %*% t(M) %*% V       # resultado: p x p if V is p x k? -> see check above
  
  # construir M: 1 en diagonal, 2 en triángulo inferior, 0 en triángulo superior
  p <- nrow(A)
  M <- matrix(0, nrow = p, ncol = p)
  diag(M) <- 1
  M[lower.tri(M)] <- 2
  
  # multiplicación elemento a elemento
  B <- A * M
  
  if (return_matrix) {
    return(B)   # si se quiere la matriz completa
  }
  
  # vech: tomar triángulo inferior por columnas (incluye diagonal)
  v <- B[lower.tri(B, diag = TRUE)]
  return(as.numeric(v))
}


pca_coef_Upsilon1_normalizado <- function( fd_centered,
                                           y,
                                           beta,
                                           gamma,
                                           alpha,
                                           step_gradient,
                                           iterations,
                                           K = NULL,                    # <<< NEW
                                           var_threshold = 0.95,
                                           basis,
                                           LdPenalization,
                                           lambda_lin = 1e-12, 
                                           lambda_quad = 1e-12, 
                                           modelo_quad = TRUE,
                                           tol = 1e-6,
                                           batch_size = 100,
                                           scale = TRUE,
                                           verbose = TRUE) 
{   
  
  cost_history <- numeric(iterations)
  psi <- inprod(basis, basis)
  
  M <- chol(psi)
  Minv <- solve(t(M))
  
  psi <- Minv %*% psi %*% t(Minv)
  
  A <- t(fd_centered$coefs)
  A_psi <- A %*% t(M) %*% psi
  
  pca_result <- prcomp(A_psi, center=TRUE, scale.=scale)
  
  ## seleccionar K si no se pasa
  if (is.null(K)) {
    var_exp_acum <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
    K <- which(var_exp_acum >= var_threshold)[1]
    
    if (is.na(K)) {
      K <- ncol(pca_result$rotation)
      warning(sprintf("No se alcanzó %.2f de varianza — usando K=%d",
                      var_threshold, K))
    }
  }
  
  if (K <= 0) stop("K debe ser >= 1")
  if (K > ncol(pca_result$rotation)) stop("K excede cantidad de PCs disponibles")
  
  Z_lin <- pca_result$x[, 1:K, drop=FALSE]
  Z_quad <- generar_matriz_cuadratica(Z_lin)
  
  LdBasis <- getbasispenalty(basis, LdPenalization)
  
  V <- pca_result$rotation[, 1:K, drop=FALSE]
  
  penalty_beta_reduced <- LdBasis / norm(LdBasis, "F")
  penalty_gamma_reduced <- LdBasis / norm(LdBasis, "F")
  
  res_gradiente <- gradient_descent_penalized_quad(
    Z_lin = Z_lin,
    Z_quad = Z_quad,
    y = y,
    transformada_beta_fun = transformada_beta_pca_normalizado,
    transformada_gamma_fun = transformada_gamma_pca_normalizado,
    H_beta_fun = H_beta_1_normalizado,
    H_gamma_fun = H_gamma_1_normalizado,
    modelo_quad = modelo_quad,
    step_gradient = step_gradient,
    tol = tol,
    batch_size = batch_size,
    lambda_lin = lambda_lin,
    lambda_quad = lambda_quad,
    iterations = iterations,
    penalty_beta_reduced = penalty_beta_reduced,
    penalty_gamma_reduced = penalty_gamma_reduced,
    beta = beta,
    gamma = gamma,
    V_lin = V,
    V_quad = V,
    psi = M,
    verbose = verbose)
  
  mu <- pca_result$center
  s <- if (scale) pca_result$scale else rep(1,length(mu))
  s_safe <- ifelse(s==0,1,s)
  
  V <- (pca_result$rotation * 1/s_safe)[,1:K,drop=FALSE]
  
  gamma_original <- H_gamma_1_normalizado(res_gradiente$UpsGamma, V, M)
  beta_original  <- H_beta_1_normalizado(res_gradiente$UpsBeta, V, M)
  
  alpha_original <- res_gradiente$alpha
  
  return(list(
    beta_reduced = res_gradiente$UpsBeta,
    gamma_reduced = res_gradiente$UpsGamma,
    beta = as.vector(beta_original),
    gamma = gamma_original,
    alpha = alpha_original,
    cost_history = res_gradiente$cost_history,
    K = K,
    pca_model = pca_result
  ))
}


############################ FPLSR #####################
##########Funciones auxiliares####################

#lo uso al hacer las plsr
logit_regression <- function(Y, X, cuad = TRUE) {
  X <- as.data.frame(X)

  model <- glm(Y ~ ., data = X, family = binomial(link = "logit"))

  coef_summary <- summary(model)$coefficients
  slope <- coef_summary["X", "Estimate"]
  se <- coef_summary["X", "Std. Error"]
  return(list(slope = slope, se = se))
}


fplslr <- function(Y, H_mat, alpha = 0.05, max_components = 5, scale = TRUE, cuad = TRUE) {
  n <- length(Y)
  p <- ncol(H_mat)
  
  # --- Normalización ---
  if (scale) {
    means <- colMeans(H_mat, na.rm = TRUE)
    sds   <- apply(H_mat, 2, sd, na.rm = TRUE)
    sds[sds == 0] <- 1  # para evitar división por cero
    
    H_mat_scaled <- sweep(H_mat, 2, means, "-")
    H_mat_scaled <- sweep(H_mat_scaled, 2, sds, "/")
  } else {
    means <- rep(0, p)
    sds   <- rep(1, p)
    H_mat_scaled <- H_mat
  }
  
  T_mat <- matrix(0, nrow = n, ncol = max_components)
  V_mat <- matrix(0, nrow = p, ncol = max_components)
  
  for (l in 1:max_components) {
    delta <- numeric(p)
    se_delta <- numeric(p)
    
    # Paso 1: para cada columna Hj, ajustar regresión logística univariada
    for (j in 1:p) {
      Hj <- H_mat_scaled[, j]
      if (l > 1) {
        Hj <- residuals(lm(Hj ~ T_mat[, 1:(l - 1)]))
      }
      reg <- logit_regression(Y, Hj, cuad = cuad)
      delta[j] <- reg$slope
      se_delta[j] <- reg$se
    }
    
    # Paso 2: construir vector de pesos v^(l)
    z_crit <- qnorm(1 - alpha / 2)
    v <- delta
    v[abs(delta / se_delta) <= z_crit] <- 0
    if (all(v == 0)) {
      message("Ningún coeficiente significativo, se detiene el algoritmo.")
      break
    }
    v <- v / sqrt(sum(v^2))  # normalizar
    V_mat[, l] <- v
    
    # Paso 3: construir componente T_l
    R_mat <- H_mat_scaled
    if (l > 1) {
      R_mat <- apply(H_mat_scaled, 2, function(hj) residuals(lm(hj ~ T_mat[, 1:(l - 1)])))
    }
    T_l <- R_mat %*% v
    T_mat[, l] <- T_l
  }
  
  # Eliminar columnas nulas de V_mat
  cols_to_keep <- colSums(V_mat^2) > 0
  V_mat <- V_mat[, cols_to_keep, drop = FALSE]
  
  # Recortar T_mat a la misma cantidad de columnas que V_mat
  T_mat <- T_mat[, cols_to_keep, drop = FALSE]
  
  return(list(
    tt = T_mat,
    ww = V_mat,
    center = means,   # <--- para des-normalizar
    scale = sds        # <--- para des-normalizar
  ))
}

################################################# Upsilon 2 ###########################################


# ---- Operador H_beta^{(2)}(Upsilon) = V_LIN^{(2)} * Upsilon ----
H_beta_2 <- function(Upsilon, V_lin) {
  V_lin <- as.matrix(V_lin)
  Upsilon <- as.vector(Upsilon)
  
  return(V_lin %*% Upsilon)
}

# ---- Operador H_gamma^{(2)}(Upsilon) =  V_lin^{(2)} %*% vech^{-1}(Upsilon * M^{-1}) ----
H_gamma_2 <- function(Upsilon, V_quad) {
  
  V_quad <- as.matrix(V_quad)
  Upsilon <- as.vector(Upsilon)
  
  # tamaño implícito de matriz simétrica
  p <- ncol(V_quad)
  len_expected <- p * (p + 1) / 2
  
  if (length(Upsilon) != len_expected) {
    stop(sprintf("Longitud de Upsilon incorrecta: esperada %d (para p=%d).", len_expected, p))
  }
  
  # construir (M/2) 
  M <- M_matrix_inv(p)
  M_vec <- M[lower.tri(M, diag = TRUE)] 
  
  # reconstruir matriz simétrica
  Gamma <- vech_inv(Upsilon * M_vec, p)
  
  # aplicar la transformación
  return(V_quad %*% Gamma %*% t(V_quad))
  
}


# vT_beta: calcula t(V) %*% beta
transformada_beta_plsr2 <- function(beta, V_lin) {
  V_lin <- as.matrix(V_lin)
  beta <- as.matrix(beta)
  
  # chequeo de compatibilidad de dimensiones
  if (nrow(V_lin) != nrow(beta)) {
    stop(sprintf("Dimensiones incompatibles: nrow(V) = %d, nrow(beta) = %d", nrow(V_lin), nrow(beta)))
  }
  
  # producto V^T * beta
  result <- t(V_lin) %*% beta
  return(result)
}



# vech_transform: devuelve vech( (t(V) %*% gamma %*% V)  *  M )
transformada_gamma_plsr2 <- function(gamma, V_quad, return_matrix = FALSE) {
  # coerciones y checks básicos
  V_quad <- as.matrix(V_quad)
  gamma <- as.matrix(gamma)
  if (!is.numeric(V_quad) || !is.numeric(gamma)) {
    stop("V y gamma deben ser matrices numéricas.")
  }
  # dimensiones compatibles: gamma debe ser p x p, V debe ser p x p' o p x k tal que t(V) %*% gamma %*% V tenga sentido.
  if (nrow(V_quad) != nrow(gamma)) {
    stop(sprintf("Dimensiones incompatibles: ncol(V) = %d, nrow(gamma) = %d", ncol(V), nrow(gamma)))
  }
  
  # calcular A = t(V) %*% gamma %*% V
  A <- t(V_quad) %*% gamma %*% V_quad       # resultado: p x p if V is p x k? -> see check above
  
  # construir M: 1 en diagonal, 2 en triángulo inferior, 0 en triángulo superior
  p <- nrow(A)
  M <- matrix(0, nrow = p, ncol = p)
  diag(M) <- 1
  M[lower.tri(M)] <- 2
  
  # multiplicación elemento a elemento
  B <- A * M
  
  if (return_matrix) {
    return(B)   # si se quiere la matriz completa
  }
  
  # vech: tomar triángulo inferior por columnas (incluye diagonal)
  v <- B[lower.tri(B, diag = TRUE)]
  return(as.numeric(v))
}


PLSR_Upsilon_2 <- function(fd_centered, y, modelo_quad = TRUE,
                    beta = NULL,
                    gamma = NULL,
                    alpha = NULL,
                    step_gradient = 0.01, iterations = 1000, batch_size = 100,
                    basis, LdPenalization,
                    nt = length(basis$nbasis),
                    lambda_lin = 1e-12, lambda_quad = 1e-12,
                    nivel_significancia = 0.05,
                    scale = TRUE, tol = 1e-6,
                    verbose = TRUE) {
  # validaciones
  if (!is.numeric(y)) y <- as.numeric(y)
  n <- length(y)
  
  psi <- inprod(basis, basis)
  A <- t(fd_centered$coefs)
  A_psi <- A %*% psi
  
  # ajustar PLS lineal 
  pls_model_lin <- fplslr(Y = y,
                          H_mat = A_psi,
                          max_components = nt,
                          scale = scale,
                          alpha = nivel_significancia)
  
  # obtener Z y V y forzar formatos
  Z_lin <- as.matrix(pls_model_lin$tt)
  storage.mode(Z_lin) <- "numeric"
  if (is.null(ncol(Z_lin))) Z_lin <- matrix(Z_lin, ncol = 1)
  if (nrow(Z_lin) != n) stop("nrow(Z_lin) distinto de length(y) -> revisar fplslr output.")
  
  # asegurar que V (ww) existe y tiene columnas >= 1
  if (!is.null(pls_model_lin$ww)) {
    V_full_lin <- as.matrix(pls_model_lin$ww)
    storage.mode(V_full_lin) <- "numeric"
  } else {
    V_full_lin <- diag(ncol(Z_lin))
  }
  
  # sincronizar K con lo que realmente tenemos en Z y en ww
  K_z_lin <- ncol(Z_lin)
  K_v_lin <- ncol(V_full_lin)
  K_lin <- min(K_z_lin, K_v_lin)   # número efectivo de componentes compartidas
  if (K_lin < 1) stop("PLS no devolvió componentes válidas.")
  
  # recortar Z y V al K efectivo (mantener forma)
  Z_lin <- Z_lin[, 1:K_lin, drop = FALSE]
  V_lin <- V_full_lin[, 1:K_lin, drop = FALSE]

  # recortar Z y V al K efectivo (mantener forma)
  Z_quad <- generar_matriz_cuadratica(Z_lin)
  V_quad <- V_lin
  
  
  # rugosidad y verificar dimensiones
  LdBasis <- getbasispenalty(basis, LdPenalization)
  LdBasis_norm <- LdBasis / norm(LdBasis, type = "F")
  penalty_beta_reduced <- LdBasis_norm
  penalty_gamma_reduced <- LdBasis_norm
  
  # llamar al optimizador robusto (le pasamos Z ya recortado)
  res_gradiente <- gradient_descent_penalized_quad(
    Z_lin = Z_lin,
    Z_quad = Z_quad,
    y = y,
    transformada_beta_fun = transformada_beta_plsr2,
    transformada_gamma_fun = transformada_gamma_plsr2,
    H_beta_fun = H_beta_2,
    H_gamma_fun = H_gamma_2,
    modelo_quad = modelo_quad,
    step_gradient = step_gradient,
    iterations = iterations,
    batch_size = batch_size,
    tol = tol,
    lambda_lin = lambda_lin,
    lambda_quad = lambda_quad,
    penalty_beta_reduced = penalty_beta_reduced,
    penalty_gamma_reduced = penalty_gamma_reduced,
    beta = beta,
    gamma = gamma,
    V_lin = V_lin,
    V_quad = V_quad,
    verbose = verbose
  )
  
  # deshacer escala: cuidar s_safe y dimensiones de V
  mu_lin <- pls_model_lin$center
  s_lin  <- if (scale) pls_model_lin$scale else rep(1, length(mu_lin))
  s_safe_lin <- ifelse(s_lin == 0, 1, s_lin)
  # V puede tener p filas y K columnas; ya recortado a K
  V_scaled_lin <- (V_lin * (1 / s_safe_lin[1:nrow(V_lin)]))
  V_scaled_lin <- as.matrix(V_scaled_lin)
  
  
  # reconstruir en espacio original
  gamma_original  <- H_gamma_2(res_gradiente$UpsGamma ,V_scaled_lin) 
  beta_original   <- H_beta_2(res_gradiente$UpsBeta ,V_scaled_lin)
  alpha_original  <- res_gradiente$alpha
  
  return(list(
    beta_reduced = res_gradiente$UpsBeta,
    gamma_reduced = res_gradiente$UpsGamma,
    beta = as.vector(beta_original),
    gamma = gamma_original,
    alpha = alpha_original,
    cost_history = res_gradiente$cost_history,
    K_lin = K_lin,
    pls_model_lin = pls_model_lin
  ))
}


################################################## Upsilon 3 ########################################

# ---- Operador H_beta^{(3)}(Upsilon) = V_LIN^{(3)} * Upsilon ----
H_beta_3 <- function(Upsilon, V_lin) {
  V_lin <- as.matrix(V_lin)
  Upsilon <- as.vector(Upsilon)
  
  return(V_lin %*% Upsilon)
}

# ---- Operador H_gamma^{(3)}(Upsilon) =  * vech^{-1}(V^{(3)}Upsilon * M^{-1}) ----
H_gamma_3 <- function(Upsilon, V_quad) {
  V_quad <- as.matrix(V_quad)
  Upsilon <- as.vector(Upsilon)
  
  
  return(reconstruir_matriz_simetrica_div2(as.vector(V_quad %*% Upsilon)))
  
}


# vT_beta: calcula t(V) %*% beta
transformada_beta_plsr3 <- function(beta, V_lin) {
  V_lin <- as.matrix(V_lin)
  beta <- as.matrix(beta)
  
  # chequeo de compatibilidad de dimensiones
  if (nrow(V_lin) != nrow(beta)) {
    stop(sprintf("Dimensiones incompatibles: nrow(V) = %d, nrow(beta) = %d", nrow(V_lin), nrow(beta)))
  }
  
  # producto V^T * beta
  result <- t(V_lin) %*% beta
  return(result)
}



# vech_transform: devuelve vech( (t(V) %*% gamma %*% V)  *  M )
transformada_gamma_plsr3 <- function(gamma, V_quad, return_matrix = FALSE) {
  # coerciones y checks básicos
  V_quad <- as.matrix(V_quad)
  gamma <- as.matrix(gamma)
  if (!is.numeric(V_quad) || !is.numeric(gamma)) {
    stop("V y gamma deben ser matrices numéricas.")
  }

  # construir M: 1 en diagonal, 2 en triángulo inferior, 0 en triángulo superior
  p <- nrow(gamma)
  M <- matrix(0, nrow = p, ncol = p)
  diag(M) <- 1
  M[lower.tri(M)] <- 2
  
  # multiplicación elemento a elemento
  B <- gamma * M
  
  # vech: tomar triángulo inferior por columnas (incluye diagonal)
  vechGamma <- as.numeric(B[lower.tri(B, diag = TRUE)])
  
  # calcular A = t(V) %*% vech(gamma )
  A <- t(V_quad) %*% vechGamma       
  
  
  return(as.numeric(A))
}



FPLSR_Upsilon3 <- function(fd_centered, y, modelo_quad = TRUE,
                                            beta = NULL,
                                            gamma = NULL,
                                            alpha = NULL,
                                            step_gradient = 0.01, iterations = 1000, batch_size = 100,
                                            basis, LdPenalization,
                                            nt_lin = length(basis$nbasis),
                                            nt_quad = length(basis$nbasis)*(length(basis$nbasis)+1)/2,
                                            lambda_lin = 1e-12, lambda_quad = 1e-12,
                                            nivel_significancia_lin = 0.05, nivel_significancia_quad= 0.05,
                                            scale = TRUE, tol = 1e-6,
                                            verbose = TRUE) {
  # validaciones
  if (!is.numeric(y)) y <- as.numeric(y)
  n <- length(y)
  
  psi <- inprod(basis, basis)
  A <- t(fd_centered$coefs)
  A_psi <- A %*% psi
  
  # ajustar PLS lineal 
  pls_model_lin <- fplslr(Y = y,
                      H_mat = A_psi,
                      max_components = nt_lin,
                      scale = scale,
                      alpha = nivel_significancia_lin)
  
  # obtener Z y V y forzar formatos
  Z_lin <- as.matrix(pls_model_lin$tt)
  storage.mode(Z_lin) <- "numeric"
  if (is.null(ncol(Z_lin))) Z_lin <- matrix(Z_lin, ncol = 1)
  if (nrow(Z_lin) != n) stop("nrow(Z_lin) distinto de length(y) -> revisar fplslr output.")
  
  # asegurar que V (ww) existe y tiene columnas >= 1
  if (!is.null(pls_model_lin$ww)) {
    V_full_lin <- as.matrix(pls_model_lin$ww)
    storage.mode(V_full_lin) <- "numeric"
  } else {
    V_full_lin <- diag(ncol(Z_lin))
  }
  
  # sincronizar K con lo que realmente tenemos en Z y en ww
  K_z_lin <- ncol(Z_lin)
  K_v_lin <- ncol(V_full_lin)
  K_lin <- min(K_z_lin, K_v_lin)   # número efectivo de componentes compartidas
  if (K_lin < 1) stop("PLS no devolvió componentes válidas.")
  
  # recortar Z y V al K efectivo (mantener forma)
  Z_lin <- Z_lin[, 1:K_lin, drop = FALSE]
  V_lin <- V_full_lin[, 1:K_lin, drop = FALSE]
  
  A_psi_quad <- as.matrix(generar_matriz_cuadratica(A_psi))
  # ajustar PLS CUARATICO (tu función)
  pls_model_quad <- fplslr(Y = y,
                          H_mat = A_psi_quad,
                          max_components = nt_quad,
                          scale = scale,
                          alpha = nivel_significancia_quad)
  
  # obtener Z y V y forzar formatos
  Z_quad <- as.matrix(pls_model_quad$tt)
  storage.mode(Z_quad) <- "numeric"
  if (is.null(ncol(Z_quad))) Z_quad <- matrix(Z_quad, ncol = 1)
  if (nrow(Z_quad) != n) stop("nrow(Z_quad) distinto de length(y) -> revisar fplslr output.")
  
  # asegurar que V (ww) existe y tiene columnas >= 1
  if (!is.null(pls_model_quad$ww)) {
    V_full_quad <- as.matrix(pls_model_quad$ww)
    storage.mode(V_full_quad) <- "numeric"
  } else {
    V_full_quad <- diag(ncol(Z_quad))
  }
  
  # sincronizar K con lo que realmente tenemos en Z y en ww
  K_z_quad <- ncol(Z_quad)
  K_v_quad <- ncol(V_full_quad)
  K_quad <- min(K_z_quad, K_v_quad)   # número efectivo de componentes compartidas
  if (K_quad < 1) stop("PLS no devolvió componentes válidas.")
  
  # recortar Z y V al K efectivo (mantener forma)
  Z_quad <- Z_quad[, 1:K_quad, drop = FALSE]
  V_quad <- V_full_quad[, 1:K_quad, drop = FALSE]
  

  # rugosidad y verificar dimensiones
  LdBasis <- getbasispenalty(basis, LdPenalization)
  LdBasis_norm <- LdBasis / norm(LdBasis, type = "F")
  penalty_beta_reduced <- LdBasis_norm
  penalty_gamma_reduced <- LdBasis_norm
  
  # llamar al optimizador robusto (le pasamos Z ya recortado)
  res_gradiente <- gradient_descent_penalized_quad(
    Z_lin = Z_lin,
    Z_quad = Z_quad,
    y = y,
    transformada_beta_fun = transformada_beta_plsr3,
    transformada_gamma_fun = transformada_gamma_plsr3,
    H_beta_fun = H_beta_3,
    H_gamma_fun = H_gamma_3,
    modelo_quad = modelo_quad,
    step_gradient = step_gradient,
    iterations = iterations,
    batch_size = batch_size,
    tol = tol,
    lambda_lin = lambda_lin,
    lambda_quad = lambda_quad,
    penalty_beta_reduced = penalty_beta_reduced,
    penalty_gamma_reduced = penalty_gamma_reduced,
    beta = beta,
    gamma = gamma,
    V_lin = V_lin,
    V_quad = V_quad,
    verbose = verbose
  )
  
  # deshacer escala: cuidar s_safe y dimensiones de V
  mu_lin <- pls_model_lin$center
  s_lin  <- if (scale) pls_model_lin$scale else rep(1, length(mu_lin))
  s_safe_lin <- ifelse(s_lin == 0, 1, s_lin)
  # V puede tener p filas y K columnas; ya recortado a K
  V_scaled_lin <- (V_lin * (1 / s_safe_lin[1:nrow(V_lin)]))
  V_scaled_lin <- as.matrix(V_scaled_lin)
  
  # deshacer escala: cuidar s_safe y dimensiones de V
  mu_quad <- pls_model_quad$center
  s_quad  <- if (scale) pls_model_quad$scale else rep(1, length(mu_quad))
  s_safe_quad <- ifelse(s_quad == 0, 1, s_quad)
  # V puede tener p filas y K columnas; ya recortado a K
  V_scaled_quad <- (V_quad * (1 / s_safe_quad[1:nrow(V_quad)]))
  V_scaled_quad <- as.matrix(V_scaled_quad)
  
  # reconstruir en espacio original
  gamma_original  <- H_gamma_3(res_gradiente$UpsGamma ,V_scaled_quad) 
  beta_original   <- H_beta_3(res_gradiente$UpsBeta ,V_scaled_lin)
  alpha_original  <- res_gradiente$alpha
  
  return(list(
    beta_reduced = res_gradiente$UpsBeta,
    gamma_reduced = res_gradiente$UpsGamma,
    beta = as.vector(beta_original),
    gamma = gamma_original,
    alpha = alpha_original,
    cost_history = res_gradiente$cost_history,
    K_lin = K_lin,
    K_quad = K_quad,
    pls_model_lin = pls_model_lin,
    pls_model_quad = pls_model_quad
    
  ))
}


##grid search cv ups3 

library(pROC)
random_search_cv_FPLSR_significance <- function(
    fd_centered, y, 
    basis, LdPenalization,
    modelo_quad = TRUE,
    
    # Nuevos candidatos
    nivel_significancia_lin_candidates,
    nivel_significancia_quad_candidates,
    size_step_candidates,
    
    # Fijos (pueden ser tuneables en otra capa)
    nt_lin = length(basis$nbasis),
    nt_quad = length(basis$nbasis)*(length(basis$nbasis)+1)/2,
    lambda_lin = 1e-12,
    lambda_quad = 1e-12,
    
    # Random search / CV
    n_iter = 20,
    K = 5,
    seed = 123,
    verbose = TRUE
){
  set.seed(seed)
  n <- length(y)
  
  #===== Generar combinaciones random =====#
  combos <- data.frame(
    nivel_lin  = sample(nivel_significancia_lin_candidates,  n_iter, replace = TRUE),
    nivel_quad = sample(nivel_significancia_quad_candidates, n_iter, replace = TRUE),
    size_step  = sample(size_step_candidates, n_iter, replace = TRUE)
  )
  
  # resultado final
  results <- data.frame()
  
  #===== Particionar folds =====#
  folds <- sample(rep(1:K, length.out = n), replace =  TRUE)
  
  #===== Loop principal =====#
  for(i in 1:nrow(combos)){
    
    nivel_lin_i  <- combos$nivel_lin[i]
    nivel_quad_i <- combos$nivel_quad[i]
    size_step_i  <- combos$size_step[i]
    
    if(verbose){
      cat("\n### Evaluando combinación", i, "de", nrow(combos), "###\n")
      cat("nivel_significancia_lin =", nivel_lin_i, 
          "| nivel_significancia_quad =", nivel_quad_i,
          "| size_step =", size_step_i, "\n")
    }
    
    auc_folds <- c()
    
    #====== CROSS VALIDATION ======#
    for(k in 1:K){
      test_idx  <- which(folds == k)
      train_idx <- setdiff(1:n, test_idx)
      
      fd_train <- fd(coef = fd_centered$coefs[, train_idx, drop = FALSE],
                     basis = basis)
      fd_test  <- fd(coef = fd_centered$coefs[, test_idx, drop = FALSE],
                     basis = basis)
      
      y_train <- y[train_idx]
      y_test  <- y[test_idx]
      
      #---- Ajustar el modelo con esta combinación ----#
      modelo_fit <- FPLSR_Upsilon3(
        fd_centered = fd_train,
        y = y_train,
        modelo_quad = modelo_quad,
        basis = basis,
        LdPenalization = LdPenalization,
        nt_lin = nt_lin,
        nt_quad = nt_quad,
        batch_size = 15,
        step_gradient = size_step_i,
        lambda_lin = lambda_lin,
        lambda_quad = lambda_quad,
        nivel_significancia_lin = nivel_lin_i,
        nivel_significancia_quad = nivel_quad_i,
        verbose = FALSE
      )
      
      #===== Predicción de probabilidades =====#
      y_prob_valid <- predecir_probabilidades_validacion(
        modelo_fit,
        fd_test,
        estimationBasis = basis,
        modo = "PLSR"
      )
      
      #===== Cálculo AUC–PR =====#
      pr_obj   <- get_pr(y_test, y_prob_valid)
      auc_pr_k <- pr_obj$auc.integral
      
      auc_folds[k] <- auc_pr_k
    }
    
    mean_auc <- mean(auc_folds)
    
    results <- rbind(
      results,
      data.frame(
        nivel_significancia_lin  = nivel_lin_i,
        nivel_significancia_quad = nivel_quad_i,
        size_step = size_step_i,
        mean_auc = mean_auc
      )
    )
  }
  
  #===== Ordenar por AUC =====#
  results <- results[order(-results$mean_auc), ]
  
  return(list(
    results = results,
    best_params = results[1, ]
  ))
}
