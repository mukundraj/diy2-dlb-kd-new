#!/bin/sh
oname="$(dirname $1)/config.nek.$2.xml"
#echo $oname
echo "$1"
#echo "$1" | sed 's/_test.rb$/_spec.rb/' > $oname
sed "s/512/$2/" $1 > $oname 