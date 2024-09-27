file.edit('.gitignore')
file.edit('.github/workflows/pkgdown.yaml')
file.edit('.github/workflows/release.yml')
file.edit('_pkgdown.yml')
file.edit('NAMESPACE')
file.edit('DESCRIPTION')
file.edit('README.md')
#-----
# pkgdown::build_favicons() # run once when you have your man/figures/logo.png
library(pkgdown)
library(roxygen2)
roxygenise(clean = TRUE)
#build_home()
#build_reference()
build_site()
preview_site()

# git commit and push  ------
# commit all changes
system('git add .')
system('git commit -m "version 1.0.2"')
# Create the tag
system('git tag -a v1.0.2 -m "Release version 1.0.2"')
# Push both the commit and the tag to the remote repository
system('git push')
system('git push origin v1.0.2')



#  usethis -----
library(usethis)
use_pipe()
use_description(fields = list(Language = "es"))
edit_r_profile(scope = c("user", "project"))
use_mit_license()       # need a LICENSE file
use_roxygen_md()        # use {roxygen2} for documentation and configuration
use_package_doc()       # setup a package-level manual page
use_testthat()          # setup testing infrastructure
use_test("placeholder") # setup a placeholder test file
devtools::document()             # Let {roxygen2} create NAMESPACE entries, build manual pages (and, more later on)
devtools::check()                # looking for the three "0's" that tell us we're ready to roll!
use_git()               # put the directory under git version control
git_vaccinate()         # Prevent leaking credentials and other unnecessary filesystem cruft


#  ----  buid pdf manual
library(roxygen2)
devtools::build_manual()
