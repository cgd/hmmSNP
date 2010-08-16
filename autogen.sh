#!/bin/sh

echo ""
echo "configuring hmmSNP build system"
echo "" 


mkdir -p buildutils

echo "calling autotools"
echo ""

#aclocal || exit
#autoheader  || exit
#autoconf -f || exit
#automake  -f -a -c --foreign || exit

autoreconf -f -v -i || exit


echo ""

echo "build system configuration complete, ready for:"
echo "   ./configure:"
echo "   make"

