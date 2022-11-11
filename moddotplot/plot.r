#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
require(scales)
library(RColorBrewer)
library(data.table)
library(cowplot)
library(glue)
library(parallel)
library(argparse)
require("argparse")

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ncolors <- 11 # max value is 11
########################################################################################################
make_scale <- function(vals) {
  comma(vals / 1e6)
}

make_k <- function(vals) {
  comma(vals / 1e3)
}

get_colors <- function(sdf) {
  bot <- floor(min(sdf$perID_by_events))
  top <- 100
  breaks <- unique(c(quantile(sdf$perID_by_events, probs = seq(0, 1, by = 1 / ncolors))))
  labels <- seq(length(breaks) - 1)
  # corner case of only one %id value
  if (length(breaks) == 1) {
    return(factor(rep(1, length(sdf$perID_by_events))))
  }
  return(cut(sdf$perID_by_events, breaks = breaks, labels = labels, include.lowest = TRUE))
}

read_bedpe <- function(all.files) {
  l <- lapply(all.files, fread, sep = "\t")
  df <- rbindlist(l)

  df$discrete <- get_colors(df)
  if ("#query_name" %in% colnames(df)) {
    df$q <- df$`#query_name`
    df$q_st <- df$query_start
    df$q_en <- df$query_end
    df$r <- df$reference_name
    df$r_st <- df$reference_start
    df$r_en <- df$reference_end
  }
  window <- max(df$q_en - df$q_st)
  df$first_pos <- df$q_st / window
  df$second_pos <- df$r_st / window
  return(df)
}

# get the lowest 0.1% of the data so we can not plot it
make_hist <- function(sdf) {
  scale_x_fun <- function(x) sprintf("%.1f", x)

  bot <- quantile(sdf$perID_by_events, probs = 0.001)[[1]]
  count <- nrow(sdf)
  extra <- ""
  my_scale <- comma
  if (count > 1e5) {
    extra <- "\n(thousands)"
    my_scale <- make_k
  }
  p <- ggplot(data = sdf, aes(perID_by_events, fill = discrete)) +
    geom_histogram(bins = 300) +
    theme_cowplot() +
    scale_fill_brewer(palette = "Spectral", direction = -1) +
    theme(legend.position = "none") +
    scale_y_continuous(labels = my_scale) +
    scale_x_continuous(labels = scale_x_fun) +
    coord_cartesian(xlim = c(bot, 100)) +
    xlab("% identity estimate") +
    ylab(glue("# of estimates{extra}"))
  p
}

diamond <- function(row) {
  # side_length, x, y, color) {
  # print(row)
  side_length <- as.numeric(row["window"])
  x <- as.numeric(row["w"])
  y <- as.numeric(row["z"])

  base <- matrix(c(1, 0, 0, 1, -1, 0, 0, -1), nrow = 2) * sqrt(2) / 2
  trans <- (base * side_length) + c(x, y)
  df <- as.data.frame(t(trans))
  colnames(df) <- c("w", "z")
  df$discrete <- as.numeric(row["discrete"])
  df$group <- as.numeric(row["group"])
  df
}

make_tri <- function(sdf, rname = "") {
  sdf$w <- (sdf$first_pos + sdf$second_pos)
  sdf$z <- -sdf$first_pos + sdf$second_pos
  window <- max(sdf$q_en - sdf$q_st)
  tri_scale <- max(sdf$q_st) / max(sdf$w)
  sdf$window <- max(sdf$q_en - sdf$q_st) / tri_scale
  sdf$group <- seq(nrow(sdf))
  df.d <- rbindlist(apply(sdf, 1, diamond))

  ggplot(df.d) +
    geom_polygon(aes(x = w * tri_scale, y = z * window, group = group, fill = factor(discrete))) +
    theme_cowplot() +
    scale_fill_brewer(palette = "Spectral", direction = -1) +
    scale_x_continuous(labels = make_scale, limits = c(0, NA)) +
    scale_y_continuous(labels = make_scale, limits = c(0, NA)) +
    xlab("Genomic position (Mbp)") +
    ylab("") +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    ) +
    ggtitle(rname)
}

make_dot <- function(sdf, rname = "") {
  max <- max(sdf$q_en, sdf$r_en)
  window <- max(sdf$query_end - sdf$query_start)
  ggplot(sdf) +
    geom_tile(aes(x = q_st, y = r_st, fill = discrete, height = window, width = window)) +
    theme_cowplot() +
    scale_fill_brewer(palette = "Spectral", direction = -1) +
    theme(legend.position = "none") +
    scale_x_continuous(labels = make_scale, limits = c(0, max)) +
    scale_y_continuous(labels = make_scale, limits = c(0, max)) +
    coord_fixed(ratio = 1) +
    facet_grid(r ~ q) +
    xlab("Genomic position (Mbp)") +
    ylab("") +
    ggtitle(rname)
}


make_plots <- function(r_name) {
  library(ggplot2)
  library(dplyr)
  require(scales)
  library(RColorBrewer)
  library(data.table)
  library(cowplot)
  library(glue)
  sdf <- copy(df[r == r_name & q == r_name])

  if (nrow(sdf) == 0) {
    return(NULL)
  }
  # set the colors
  if (!ONECOLORSCALE) {
    sdf$discrete <- get_colors(sdf)
  }

  # make the plots
  if (TRI) {
    p_lone <- make_tri(sdf, rname = r_name)
    scale <- 1 / 2.25
  } else {
    scale <- 1
    p_lone <- make_dot(sdf, rname = r_name)
  }
  p_hist <- make_hist(sdf)


  # setup the space
  dir.create(glue("{OUT}/pdfs/{r_name}/"))
  dir.create(glue("{OUT}/pngs/{r_name}/"))

  # save the plots
  if (TRI) {
    p_hist <- p_hist +
      theme(text = element_text(size = 8), axis.text = element_text(size = 8))

    # ranges for inset hist
    mmax <- max(sdf$q_en, sdf$r_en)
    build <- ggplot_build(p_lone)
    yr <- build$layout$panel_params[[1]]$y.range
    xmin <- -1 / 10 * mmax
    xmax <- mmax * 1 / 3.75
    ymin <- yr[2] * 1 / 2
    ymax <- yr[2] * 2 / 2
    print(paste(r_name, xmin, xmax, ymin, ymax))
    # combine
    plot <- p_lone + annotation_custom(
      ggplotGrob(p_hist),
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax
    )
  } else {
    plot <- cowplot::plot_grid(
      p_lone, p_hist,
      ncol = 1,
      rel_heights = c(3 * scale, 1)
    )
  }
  # if we are not making a facet only save the pdf
  if (!ONECOLORSCALE) {
    ggsave(
      plot = plot,
      file = glue("{OUT}/pdfs/{r_name}/{PRE}__{r_name}__tri.{TRI}.pdf"),
      height = 12 * scale,
      width = 9
    )

    ggsave(
      plot = plot,
      file = glue("{OUT}/pngs/{r_name}/{PRE}__{r_name}__tri.{TRI}.png"),
      height = 12 * scale,
      width = 9,
      dpi = DPI
    )
  }
  if (TRI) {
    return(plot)
  } else {
    return(p_lone)
  }
}
########################################################################################################
#
# EDIT THIS SECION FOR YOUR INPUTS
#
parser <- ArgumentParser()
parser$add_argument("-b", "--bed", help = "bedfile with alignment information")
parser$add_argument("-t", "--threads", type = "integer", help = "number of threads")
parser$add_argument("-p", "--prefix", help = "Prefix for the outputs")
args <- parser$parse_args()

PRE <- args$prefix
GLOB <- args$bed
OUT <- glue("results/{PRE}_figures")
print(PRE)
print(GLOB)

DPI <- 600
#
# STOP EDITING
#
########################################################################################################

dir.create(OUT)
dir.create(glue("{OUT}/pdfs"))
dir.create(glue("{OUT}/pngs"))
all.files <- Sys.glob(GLOB)
df <- read_bedpe(all.files) # [1:5e4]
# df=fread(GLOB)
Qs <- unique(df$q)
N <- length(Qs)
columns <- ceiling(sqrt(N + 1))
rows <- ceiling((N + 1) / columns)


vals <- c(TRUE, FALSE)
for (TRI in vals) {
  if (TRI) {
    scale <- 2 / 3
  } else {
    scale <- 1
  }
  for (ONECOLORSCALE in vals) {
    # plots = lapply(Qs, make_plots)
    cl <- parallel::makeCluster(min(args$threads, length(Qs)),
      setup_strategy = "sequential",
      outfile = ""
    )
    clusterExport(
      cl = cl,
      varlist = c(
        "TRI", "OUT", "df",
        "PRE", "DPI",
        "ONECOLORSCALE", "ncolors",
        "make_scale",
        "make_k",
        "get_colors",
        "read_bedpe",
        "make_hist",
        "diamond",
        "make_tri",
        "make_dot"
      )
    )
    plots <- parLapply(cl = cl, Qs, make_plots)
    stopCluster(cl)
    if (!TRI) {
      plots[[N + 1]] <- make_hist(df)
    }
    N <- length(plots)
    columns <- ceiling(sqrt(N))
    rows <- ceiling((N) / columns)
    p <- cowplot::plot_grid(plotlist = plots, nrow = rows, ncol = columns, labels = "auto")
    ggsave(glue("{OUT}/pdfs/{PRE}.tri.{TRI}__onecolorscale.{ONECOLORSCALE}__all.pdf"), plot = p, height = 6 * rows * scale, width = 6 * columns)
    ggsave(glue("{OUT}/pngs/{PRE}.tri.{TRI}__onecolorscale.{ONECOLORSCALE}__all.png"), plot = p, height = 6 * rows * scale, width = 6 * columns, dpi = DPI)
  }
}

#
# big plot
#
if (T) {
  facet_fig <- cowplot::plot_grid(make_hist(df), make_dot(df), rel_heights = c(1, 4), ncol = 1)
  ggsave(plot = facet_fig, file = glue("{OUT}/pdfs/{PRE}.facet.all.pdf"), height = 20, width = 16)
  ggsave(plot = facet_fig, file = glue("{OUT}/pngs/{PRE}.facet.all.png"), height = 20, width = 16, dpi = DPI)
}