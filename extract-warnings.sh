#!/bin/sh

if [ "x$1" = "x" ]; then
	echo "usage: extract-warnings.sh <logfile>"
	exit 1
fi

thisdir="`(cd \`dirname $0\`; /bin/pwd)`"
grep -i "$thisdir.*warning" "$1"
