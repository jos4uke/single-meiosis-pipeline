#!/usr/bin/env bash

## Copyright (c) 2013 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/.

################################################
##
## @script: install.sh
##
## @author: Joseph Tran (Joseph.Tran@versailles.inra.fr)
##
## @version: 0.0.1
##
## @date: 2014-01-20
##
## @description: This bash script performs the bash common functions library installation to the production server
##
###############################################


ARGV=("$@")
ARGC=$#
ERROR_LOG=/tmp/singlemeiosis_lib_install_error.log
libname="singlemeiosis"

if [[ $ARGC -ne 2 ]]; then
    echo -e "Incorrect number of required positional arguments. Should be equal to 2: \n- the first argument corresponding to the tag version of $libname library in git repository, \n- and the second and last argument to the prefix used to install the library." | tee $ERROR_LOG 2>&1
    #usage;
    exit 1
else
	LIB_VERSION=${ARGV[0]}
	PREFIX=${ARGV[1]}
    echo -e "$(date '+%Y-%m-%d %H:%M:%S') Starting the $libname library installation" | tee $ERROR_LOG 2>&1
    echo "$(date '+%Y-%m-%d %H:%M:%S') Will install library version tag: $LIB_VERSION" | tee -a $ERROR_LOG 2>&1
    echo "$(date '+%Y-%m-%d %H:%M:%S') Will install library in: $PREFIX" | tee -a $ERROR_LOG 2>&1    
fi 

# export $libname library from git repo 
LIB_VERSION_DIR=/usr/local/archives/$libname/$LIB_VERSION
mkdir -p $LIB_VERSION_DIR 2>$ERROR_LOG
rtrn=$?
if [[ rtrn -ne 0 ]]; then
    echo -e "$LIB_VERSION_DIR directory creation failed with exit error code $rtrn." | tee -a $ERROR_LOG 2>&1
    echo -e "You can get more information about $libname library installation in $ERROR_LOG." 2>&1
    exit $rtrn
else
    echo -e "$LIB_VERSION_DIR directory was created successfully." | tee -a $ERROR_LOG 2>&1
    echo -e "Will export $libname library project from git repository." | tee -a $ERROR_LOG 2>&1
fi
git archive --prefix=${libname}_$LIB_VERSION$PREFIX/ $LIB_VERSION | (cd $LIB_VERSION_DIR && tar xf -)
rtrn=$?
if [[ rtrn -ne 0 ]]; then
    echo -e "git archive failed to export libname library project with exit error code $rtrn." | tee -a $ERROR_LOG 2>&1
    echo -e "You can get more information about $libname library installation in $ERROR_LOG." 2>&1
    exit $rtrn
else
    echo -e "$LIB_VERSION_DIR directory was created successfully." | tee -a $ERROR_LOG 2>&1
    echo -e "Will install $libname library files to $PREFIX" | tee -a $ERROR_LOG 2>&1
fi

# install $libname library
#cp -R * $PREFIX/. 2>&1
cd  $LIB_VERSION_DIR/${libname}_$LIB_VERSION/ 2>$ERROR_LOG
#rm -rf usr/local/test 2>$ERROR_LOG
for f in $(find . -name \* -print 2>$ERROR_LOG); do
    echo -e "cp $f in /$f" | tee -a $ERROR_LOG 2>&1
    cp --parents $f / 2>$ERROR_LOG
done

# set privileges
#echo -ne "Setting executable privileges on mutdetect.sh script ..." | tee -a $ERROR_LOG 2>&1
#chmod 755 $PREFIX/bin/mutdetect.sh 2>$ERROR_LOG
#echo -e "OK" | tee -a $ERROR_LOG 2>&1

echo
echo "NOTE: You can get more information about libname library installation in $ERROR_LOG."
echo
echo
echo "************************************************************"
echo "*             INSTALLATION IS NOW COMPLETE                 *"
echo "************************************************************"
echo
echo
echo "The $libname library is now available in $PREFIX/share/${libname}-pipeline/lib ."
echo
