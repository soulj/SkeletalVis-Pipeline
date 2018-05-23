#
# One off installation script for DAVIDQuery library.
# Run as root.
#

if(require(DAVIDQuery)) {
  remove.packages("DAVIDQuery")
}
install.packages('DAVIDQuery.tgz', repos=NULL, type="source")
