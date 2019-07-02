#!/bin/bash

## PPR Tag for build
PPR_HASH=dbd8868
## Subdirectory in PPR repo to use
PPR_SUBDIR=src

## Available protocols for acquiring PPR source code
PPR_GIT_HTTP_ADDRESS=https://github.com/dengwirda/PPR.git
PPR_GIT_SSH_ADDRESS=git@github.com:dengwirda/PPR.git

GIT=`which git`
SVN=`which svn`
PROTOCOL=""

# PPR exists. Check to see if it is the correct version.
# Otherwise, flush the directory to ensure it's updated.
if [ -d ppr ]; then

	if [ -d .ppr_all/.git ]; then
		cd .ppr_all
		CURR_TAG=$(git describe --tags)
		cd ../
		if [ "${CURR_TAG}" == "${PPR_HASH}" ]; then
			echo "PPR version is current. Skip update"
		else
			unlink ppr
			rm -rf .ppr_all
		fi
	else
		unlink ppr
		rm -rf .ppr_all
	fi
fi


# PPR Doesn't exist, need to acquire souce code
if [ ! -d ppr ]; then
	if [ -d .ppr_all ]; then
		rm -rf .ppr_all
	fi

	if [ "${GIT}" != "" ]; then
		echo " ** Using git to acquire ppr source. ** "
		PROTOCOL="git ssh"
		git clone ${PPR_GIT_HTTP_ADDRESS} .ppr_all &> /dev/null
		if [ -d .ppr_all ]; then
			cd .ppr_all
			git checkout ${PPR_HASH} &> /dev/null
			cd ../
			ln -sf .ppr_all/${PPR_SUBDIR} ppr
		else
			git clone ${PPR_GIT_HTTP_ADDRESS} .ppr_all &> /dev/null
			PROTOCOL="git http"
			if [ -d .ppr_all ]; then
				cd .ppr_all
				git checkout ${PPR_HASH} &> /dev/null
				cd ../
				ln -sf .ppr_all/${PPR_SUBDIR} ppr
			fi
		fi
	fi
fi

if [ ! -d ppr ]; then
	echo " ****************************************************** "
	echo " ERROR: Build failed to acquire PPR source."
	echo ""
	echo " Please ensure your proxy information is setup properly for"
	echo " the protocol you use to acquire PPR."
	echo ""
	echo " The automated script attempted to use: ${PROTOCOL}"
	echo ""
	if [ "${PROTOCOL}" == "git http" ]; then
		echo " This protocol requires setting up the http.proxy git config option."
	elif [ "${PROTOCOL}" == "git ssh" ]; then
		echo " This protocol requires having ssh-keys setup, and ssh access to git@github.com."
		echo " Please use 'ssh -vT git@github.com' to debug issues with ssh keys."
	elif [ "${PROTOCOL}" == "svn" ]; then
		echo " This protocol requires having svn proxys setup properly in ~/.subversion/servers."
	elif [ "${PROTOCOL}" == "wget" ]; then
		echo " This protocol requires having the http_proxy and https_proxy environment variables"
		echo " setup properly for your shell."
	fi
	echo ""
	echo " ****************************************************** "
fi

cp Makefile-ppr ppr/Makefile
