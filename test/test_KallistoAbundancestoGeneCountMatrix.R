#test for the command line based arguments
output<-tempfile()
species<-"Mouse"
filepaths<-"../test/resources/RNASeq/testAbundance1.abundances ../test/resources/RNASeq/testAbundance2.abundances"
filenames<-"test1 test2"

command<-paste("Rscript '../alt-src/kallisto/KallistoAbundancestoGeneCountMatrix.R' --output",output,"--species",species,"--filepaths",filepaths,"--filenames",filenames,sep=" ")

system(command)

testFile<-load(output)
str(get(testFile))
