println <- function(...) cat(sprintf(...), '\n')

# Invoked by Jenkins to:
#  - define the function
#  - run unit test
#  - generate the galaxy tool

# ----------------- PARSE ARGS -----------------
# ARGS REQUIRED:
#  - functionName
#  - toolSectionName
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop('Require arguments for function and tool name.')
}
functionName <- args[1]
toolSectionName <- args[2]
println('Load tool in section "%s" using function: %s\n', toolSectionName, functionName)

# ----------------- LOAD FUNCTION -----------------
setwd(paste('src/', toolSectionName, sep=''))
println('Now in: %s', getwd())
file <- paste(functionName, '.R', sep='')
println('Loading from file: %s', file)
source(file)

# ----------------- TEST FUNCTION -----------------
println('Run function test.\n')
#TODO

# ----------------- INSTALL FUNCTION -----------------
println('Begin Galaxy tool generation.\n')
# Build in the tmp directory.
# Must copy man page so it can be found in the working directory.
setwd('../../tmp')
invisible(file.copy(paste('../man/', functionName, sep=''), '.'))

# Currently no reason why we'd need the tool ID to differ from the function name.
toolId <- functionName

if (!exists("galaxyHome")) {
  galaxyHome <- getwd()
}
galaxy(
 functionName,
 version="1.0",
 galaxyConfig=GalaxyConfig(galaxyHome, functionName, toolSectionName, toolId)
)
