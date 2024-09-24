rmarkdown::find_pandoc()


#install dependent pkgs first
to.install <- c("plotrix","ellipse","Hmisc","gplots","fields","RColorBrewer","colorspace","mnormt","Deriv","tidyr","dplyr","ggplot2","viridis", "abind", "rmarkdown", "pander", "kableExtra")
new.packages <- to.install[!(to.install %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#install WHAM
devtools::install_github("timjmiller/wham", dependencies=TRUE, ref = "devel", INSTALL_opts=c("--no-multiarch"))

#TEST
library(wham)
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
input <- prepare_wham_input(asap3) 
nofit <- fit_wham(input, do.fit = FALSE)
plot_wham_output(nofit)


#workshop git repo
git@github.com:timjmiller/wham_workshop_MUN_2024.git

usethis::use_git()
