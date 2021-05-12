
################################### 
# store_packages.R
#
# stores a list of your currently installed packages
# uncomment to generate '.rda file'

#tmp = installed.packages()

#installedpackages = as.vector(tmp[is.na(tmp[,"Priority"]), 1])
#save(installedpackages, file = 'installed_packages.rda')

####################################
# restore_packages.R
#
# installs each package from the stored list of packages

load("/home/analysis/installed_packages.rda")

for (count in 1:length(installedpackages)) install.packages(installedpackages[count])
