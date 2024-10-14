file.edit('.gitignore')
file.edit('.github/workflows/pkgdown.yaml')
file.edit('.github/workflows/Release.yml')
file.edit('.github/workflows/R-CMD-check.yml')
file.edit('_pkgdown.yml')
file.edit('NAMESPACE')
file.edit('DESCRIPTION')
file.edit('README.md')
file.edit('NEWS.Rmd')
# build website -----
# pkgdown::build_favicons() # run once when you have your man/figures/logo.png
library(devtools)
library(pkgdown)
library(roxygen2)
roxygenise(clean = TRUE)
#build_home()
#build_reference()
build_site()
preview_site()
#  load developing package  -------
devtools::load_all()


# git commit and push  ------
rmarkdown::render("NEWS.Rmd", output_file = "NEWS.md")
toolkit::git_commit_push("version 1.0.5")


#  ----  buid pdf manual
library(roxygen2)
devtools::build_manual()
