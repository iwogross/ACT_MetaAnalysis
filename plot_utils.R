# Copied and slightly altered functions from ggridges package
## Necessary to hijack plotting of quantile lines on density plots and instead
## plot 95% HPD Credible Intervals and optionally posterior modes
### (from `coda::HPDinterval()` and `MCMCglmm::posterior.mode()`)
stat_density_ridges_HPDCrI <- function(mapping = NULL, data = NULL, geom = "density_ridges",
                     position = "identity", na.rm = FALSE, show.legend = NA,
                     inherit.aes = TRUE, bandwidth = NULL, from = NULL, to = NULL,
                     jittered_points = FALSE, quantile_lines = FALSE, calc_ecdf = FALSE, quantiles = 4,
                     quantile_fun = quantile, n = 512, ...)
{
  layer(
    stat = StatDensityRidgesHPDCrI,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(bandwidth = bandwidth,
                  from = from,
                  to = to,
                  calc_ecdf = calc_ecdf,
                  quantiles = quantiles,
                  jittered_points = jittered_points,
                  quantile_lines = quantile_lines,
                  quantile_fun = quantile_fun,
                  n = n,
                  na.rm = na.rm, ...)
  )
}


################################################################################
StatDensityRidgesHPDCrI <- ggproto("StatDensityRidgesHPDCrI", Stat,
  required_aes = "x",

  default_aes = aes(height = ..density..),

  calc_panel_params = function(data, params) {
    if (is.null(params$bandwidth)) {
      xdata <- na.omit(data.frame(x=data$x, group=data$group))
      xs <- split(xdata$x, xdata$group)
      xs_mask <- vapply(xs, length, numeric(1)) > 1
      bws <- vapply(xs[xs_mask], bw.nrd0, numeric(1))
      bw <- mean(bws, na.rm = TRUE)
      message("Picking joint bandwidth of ", signif(bw, 3))

      params$bandwidth <- bw
    }

    if (is.null(params$from)) {
      params$from <- min(data$x, na.rm=TRUE) - 3 * params$bandwidth
    }

    if (is.null(params$to)) {
      params$to <- max(data$x, na.rm=TRUE) + 3 * params$bandwidth
    }

    data.frame(
      bandwidth = params$bandwidth,
      from = params$from,
      to = params$to
    )
  },

  setup_params = function(self, data, params) {
    # calculate bandwidth, min, and max for each panel separately
    panels <- split(data, data$PANEL)
    pardata <- lapply(panels, self$calc_panel_params, params)
    pardata <- ggridges:::reduce(pardata, rbind)

    if (length(params$quantiles) > 1 &&
        (max(params$quantiles, na.rm = TRUE) > 1 || min(params$quantiles, na.rm = TRUE) < 0)) {
      stop('invalid quantiles used: c(', paste0(params$quantiles, collapse = ','), ') must be within [0, 1] range')
    }

    params$bandwidth <- pardata$bandwidth
    params$from <- pardata$from
    params$to <- pardata$to
    params
  },

  compute_group = function(data, scales, from, to, bandwidth = 1,
                           calc_ecdf = FALSE, jittered_points = FALSE, quantile_lines = FALSE,
                           quantiles = 4, quantile_fun = quantile, n = 512) {
    # ignore too small groups
    if(nrow(data) < 3) return(data.frame())

    if (is.null(calc_ecdf)) calc_ecdf <- FALSE
    if (is.null(jittered_points)) jittered_points <- FALSE
    if (is.null(quantile_lines)) quantile_lines <- FALSE

    # when quantile lines are requested, we also calculate ecdf
    # this simplifies things for now; in principle, could disentangle
    # the two
    if (quantile_lines) calc_ecdf <- TRUE

    panel <- unique(data$PANEL)
    if (length(panel) > 1) {
      stop("Error: more than one panel in compute group; something's wrong.")
    }
    panel_id <- as.numeric(panel)

    d <- stats::density(
      data$x,
      bw = bandwidth[panel_id], from = from[panel_id], to = to[panel_id], na.rm = TRUE,
      n = n
    )

    # calculate maximum density for scaling
    maxdens <- max(d$y, na.rm = TRUE)

    # make interpolating function for density line
    densf <- approxfun(d$x, d$y, rule = 2)

    # calculate jittered original points if requested
    if (jittered_points) {
      df_jittered <- data.frame(
        x = data$x,
        # actual jittering is handled in the position argument
        density = densf(data$x),
        ndensity = densf(data$x) / maxdens,
        datatype = "point", stringsAsFactors = FALSE)

      # see if we need to carry over other point data
      # capture all data columns starting with "point", as those are relevant for point aesthetics
      df_points <- data[grepl("point_", names(data))]

      # uncomment following line to switch off carrying over data
      #df_points <- data.frame()

      if (ncol(df_points) == 0) {
        df_points <- NULL
        df_points_dummy <- NULL
      }
      else {
        # combine additional points data into results dataframe
        df_jittered <- cbind(df_jittered, df_points)
        # make a row of dummy data to merge with the other dataframes
        df_points_dummy <- na.omit(df_points)[1, , drop = FALSE]
      }
    } else {
      df_jittered <- NULL
      df_points_dummy <- NULL
    }

    # calculate quantiles, needed for both quantile lines and ecdf
    if ((length(quantiles)==1) && (all(quantiles >= 1))) {
      if (quantiles > 1) {
        probs <- seq(0, 1, length.out = quantiles + 1)[2:quantiles]
      }
      else {
        probs <- NA
      }
    } else {
      probs <- quantiles
      probs[probs < 0 | probs > 1] <- NA
    }

    qx <- na.omit(quantile_fun(data$x, probs = probs))

    # if requested, add data frame for quantile lines
    df_quantiles <- NULL

    if (quantile_lines && length(qx) > 0) {
      qy <- densf(qx)
      df_quantiles <- data.frame(
        x = qx,
        density = qy,
        ndensity = qy / maxdens,
        datatype = "vline",
        stringsAsFactors = FALSE
      )
      if (!is.null(df_points_dummy)){
        # add in dummy points data if necessary
        df_quantiles <- data.frame(df_quantiles, as.list(df_points_dummy))
      }
    }

    # combine the quantiles and jittered points data frames into one, the non-density frame
    df_nondens <- rbind(df_quantiles, df_jittered)

    if (calc_ecdf) {
      n <- length(d$x)
      ecdf <- c(0, cumsum(d$y[1:(n-1)]*(d$x[2:n]-d$x[1:(n-1)])))
      ecdf_fun <- approxfun(d$x, ecdf, rule = 2)
      ntile <- findInterval(d$x, qx, left.open = TRUE) + 1 # if make changes here, make them also below

      if (!is.null(df_nondens)) {
        # we add data for ecdf and quantiles back to all other data points
        df_nondens <- data.frame(
          df_nondens,
          ecdf = ecdf_fun(df_nondens$x),
          quantile = findInterval(df_nondens$x, qx, left.open = TRUE) + 1
        )
      }

      df_density <- data.frame(
        x = d$x,
        density = d$y,
        ndensity = d$y / maxdens,
        ecdf = ecdf,
        quantile = ntile,
        datatype = "ridgeline",
        stringsAsFactors = FALSE
      )
    }
    else {
      df_density <- data.frame(
        x = d$x,
        density = d$y,
        ndensity = d$y / maxdens,
        datatype = "ridgeline",
        stringsAsFactors = FALSE
      )
    }

    if (!is.null(df_points_dummy)){
      # add in dummy points data if necessary
      df_density <- data.frame(df_density, as.list(df_points_dummy))
    }

    # now combine everything and turn quantiles into factor
    df_final <- rbind(df_density, df_nondens)
    if ("quantile" %in% names(df_final)) {
      df_final$quantile <- factor(df_final$quantile)
    }

    df_final
  }
)


################################################################################
# user defined functions for argument "quantile_fun" in geom_density_ridges(2)
## 95% HPD Credible Intervals and posterior mode
HPDpstMode_fun <- function(x, probs){
  qx <- c(HPDinterval(as.mcmc(x), prob = probs)[1, , drop = TRUE],
    posterior.mode(as.mcmc(x)))
  qx <- qx[c(1,3,2)]
  names(qx)[2] <- "mode"
 qx
}    
## Just 95% HPD Credible Intervals
HPD_fun <- function(x, probs){
  HPDinterval(as.mcmc(x), prob = probs)[1, , drop = TRUE]
}

