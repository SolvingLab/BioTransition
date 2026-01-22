pkgname <- "BioTransition"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "BioTransition-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('BioTransition')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("LDNB")
### * LDNB

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LDNB
### Title: Landscape DNB Analysis
### Aliases: LDNB

### ** Examples

## No test: 
# LDNB requires a reference state
# See vignette for detailed examples
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LDNB", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LcDNB")
### * LcDNB

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LcDNB
### Title: Local Conventional DNB Analysis
### Aliases: LcDNB

### ** Examples

## No test: 
# See vignette for detailed examples
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LcDNB", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cDNB")
### * cDNB

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cDNB
### Title: Conventional DNB analysis
### Aliases: cDNB

### ** Examples

# Create example data
set.seed(42)
n_genes <- 100
n_samples <- 15

expr <- matrix(
  rnorm(n_genes * n_samples, mean = 10, sd = 2),
  nrow = n_genes, ncol = n_samples
)
rownames(expr) <- paste0("Gene", seq_len(n_genes))
colnames(expr) <- paste0("Sample", seq_len(n_samples))

state <- data.frame(
  sample_id = colnames(expr),
  state = rep(c("A", "B", "C"), each = 5)
)

## No test: 
result <- cDNB(
  expr = expr,
  state = state,
  state.levels = c("A", "B", "C"),
  min.size = 5,
  max.size = 50
)
result$DNB.score
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cDNB", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ppi_h")
### * ppi_h

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ppi_h
### Title: Human Protein-Protein Interaction Network
### Aliases: ppi_h
### Keywords: datasets

### ** Examples

# Load human PPI
data("ppi_h")

# Explore the network
dim(ppi_h)
head(ppi_h)

# Check gene coverage
all_genes <- unique(c(ppi_h$G1, ppi_h$G2))
length(all_genes)

# Filter by confidence
high_conf <- ppi_h[ppi_h$combined_score >= 700, ]
nrow(high_conf)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ppi_h", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ppi_m")
### * ppi_m

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ppi_m
### Title: Mouse Protein-Protein Interaction Network
### Aliases: ppi_m
### Keywords: datasets

### ** Examples

# Load mouse PPI
data("ppi_m")

# Explore the network
dim(ppi_m)
head(ppi_m)

# Check gene coverage
all_genes <- unique(c(ppi_m$G1, ppi_m$G2))
length(all_genes) # 21,645 genes

# Filter by confidence
high_conf <- ppi_m[ppi_m$combined_score >= 700, ]
nrow(high_conf)

## Not run: 
##D # Use in DNB analysis
##D result <- MDNB(
##D   expr = mouse_expr,
##D   state = sample_groups,
##D   state.levels = c("Control", "Treatment"),
##D   ppi = ppi_m,
##D   min.combined.score = 700
##D )
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ppi_m", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("tDNB")
### * tDNB

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: tDNB
### Title: Topological DNB analysis
### Aliases: tDNB

### ** Examples

# Create example data
set.seed(42)
n_genes <- 100
n_samples <- 15

expr <- matrix(
  rnorm(n_genes * n_samples, mean = 10, sd = 2),
  nrow = n_genes, ncol = n_samples
)
rownames(expr) <- paste0("Gene", seq_len(n_genes))
colnames(expr) <- paste0("Sample", seq_len(n_samples))

state <- data.frame(
  sample_id = colnames(expr),
  state = rep(c("A", "B", "C"), each = 5)
)

## No test: 
result <- tDNB(
  expr = expr,
  state = state,
  state.levels = c("A", "B", "C"),
  min.size = 5,
  max.size = 50
)
result$DNB.score
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("tDNB", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
