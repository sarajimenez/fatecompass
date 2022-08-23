# libraries ---------------------------------------------------------------
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(shiny)
library(shinydashboard)

# Read data mouse ---------------------------------------------------------
u2_m  <- read.table(file = "www/data/mouse/umap2_E15_5.txt")
colnames(u2_m) <- c("umap1","umap2")
i_m   <- read.table(file = "www/data/mouse/mycellID_E15_5.txt", check.names = FALSE, sep = '\t')
g_m   <- read.table(file = "www/data/mouse/g_m.txt", header = TRUE)
a_m   <- as.data.frame(read.table(file = "www/data/mouse/Areg.csv", header = TRUE, sep = ';'))
exp_traj_alpha_m <- read.table(file = "www/data/mouse/exp_traj_alpha.txt")
exp_traj_beta_m <- read.table(file = "www/data/mouse/exp_traj_beta.txt")
act_traj_alpha_m <- read.table(file = "www/data/mouse/act_traj_alpha.txt")
act_traj_beta_m <- read.table(file = "www/data/mouse/act_traj_beta.txt")

colnames(g_m) <- str_replace(colnames(g_m), "-", ".")
genes_m       <- colnames(g_m)
motifs_m      <- colnames(a_m)

dma_m <- read.table(file = "www/data/mouse/dma.csv", header = TRUE, sep = ',')
dma_m[,5] <- c()

# Parameters mouse --------------------------------------------------------

colnames(exp_traj_alpha_m)  <- colnames(g_m)
colnames(exp_traj_beta_m)   <- colnames(g_m)
colnames(act_traj_alpha_m)  <- colnames(a_m)
colnames(act_traj_beta_m)   <- colnames(a_m)

cols_m <- c("Pre-endocrine" = "#EE792F",
            "Ductal"        = "#8EBC8F",
            "Alpha"         = "#357AC5",
            "Ngn3 high EP"  = "#F5C172",
            "Delta"         = "#6A3B99",
            "Beta"          = "#A5E17F",
            "Ngn3 low EP"   = "#F1A05A",
            "Epsilon"       = "#C9B2D5")

cols_e_m <- c("Pre-endocrine" = "#c9c9c9",
              "Ductal"        = "#c9c9c9",
              "Alpha"         = "#357AC5",
              "Ngn3 high EP"  = "#c9c9c9",
              "Delta"         = "#c9c9c9",
              "Beta"          = "#A5E17F",
              "Ngn3 low EP"   = "#c9c9c9",
              "Epsilon"       = "#c9c9c9")

# Read data human ---------------------------------------------------------
u2_h  <- read.table(file = "www/data/human/umap2.txt")
colnames(u2_h) <- c("umap1","umap2")
i_h   <- read.table(file = "www/data/human/cellid.csv", header = TRUE, check.names = FALSE, sep = ',')
g_h   <- read.table(file = "www/data/human/g_h.txt", header = TRUE)
a_h   <- as.data.frame(read.table(file = "www/data/human/Areg_4.csv", header = TRUE, sep = ';'))
exp_traj_alpha_h <- read.table(file = "www/data/human/exp_traj_alpha.txt")
exp_traj_beta_h <- read.table(file = "www/data/human/exp_traj_beta.txt")
exp_traj_ec_h <- read.table(file = "www/data/human/exp_traj_ec.txt")
act_traj_alpha_h <- read.table(file = "www/data/human/act_traj_alpha.txt")
act_traj_beta_h <- read.table(file = "www/data/human/act_traj_beta.txt")
act_traj_ec_h <- read.table(file = "www/data/human/act_traj_ec.txt")

colnames(g_h) <- str_replace(colnames(g_h), "-", ".")
genes_h       <- colnames(g_h)
motifs_h      <- colnames(a_h)

dma_h <- read.table(file = "www/data/human/dma.csv", header = TRUE, sep = ',')
dma_h[,5] <- c()

# Parameters human --------------------------------------------------------

colnames(exp_traj_alpha_h)  <- colnames(g_h)
colnames(exp_traj_beta_h)   <- colnames(g_h)
colnames(exp_traj_ec_h)   <- colnames(g_h)
colnames(act_traj_alpha_h)  <- colnames(a_h)
colnames(act_traj_beta_h)   <- colnames(a_h)
colnames(act_traj_ec_h)   <- colnames(a_h)

cols_h <- c("prog_nkx61"       = "#66C3A6",
            "neurog3_early"    = "#F3C27B",
            "neurog3_mid"      = "#EE8532",
            "neurog3_late"     = "#9B5933",
            "fev_high_isl_low" = "#C9B1D7",
            "sst_hhex"         = "#6B389A",
            "Alpha-like"       = "#357AC5",
            "Beta-like"        = "#A5E17F",
            "EC-like"          = "#EE847C")

cols_e_h <- c("prog_nkx61"       = "#c9c9c9",
              "neurog3_early"    = "#c9c9c9",
              "neurog3_mid"      = "#c9c9c9",
              "neurog3_late"     = "#c9c9c9",
              "fev_high_isl_low" = "#c9c9c9",
              "sst_hhex"         = "#c9c9c9",
              "Alpha-like"       = "#357AC5",
              "Beta-like"        = "#A5E17F",
              "EC-like"          = "#EE847C")





