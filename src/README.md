Source R scripts reside in this directory.

The sub directories correspond to groups of tools within Galaxy and R scripts in these directories will get loaded into Galaxy for users to use.

To create a new R script:

*  Add a new file to an existing sub directory, or create new directory. The file must have a .R extension.
* Create the man page for the script in the man directory. This file must have te same name as the R script but with no extension.
* Commit the files to git and wait for Jenkins to checkout the project and deploy to the Galaxy server (or manually run the systems-biology-galaxy-r job in Jenkins)
