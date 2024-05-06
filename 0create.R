file.edit('.gitignore')
#file.edit('.github/workflows/_pkgdown.yml')
file.edit('_pkgdown.yml')
file.edit('NAMESPACE')
file.edit('DESCRIPTION')
#-----
library(pkgdown)
#usethis::use_github_action("pkgdown")
#usethis::use_pkgdown()
roxygen2::roxygenise(clean = TRUE)
build_home()
build_site()
preview_site()
# run once when you have your man/figures/logo.png
build_favicons(overwrite = TRUE)
#  usethis -----
library(usethis)
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
