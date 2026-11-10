## Test Environments

* Local MacOS 15.7.4 install, R 4.6
* macos-latest,   r: 'release'
* windows-latest, r: 'release'
* ubuntu-latest,   r: 'devel', http-user-agent: 'release'
* ubuntu-latest,   r: 'release'
* ubuntu-latest,   r: 'oldrel-1'

## R CMD check results

── R CMD check results ───────────────────────────────────────────────────────── profExtrema 0.2.2 ────
Duration: 2m 26.3s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Update of exisiting CRAN package

* This is a patch on an existing CRAN package to address the NOTES currently appearing in the CRAN checks. 
* The package passes devtools::check() on all test environments. 
