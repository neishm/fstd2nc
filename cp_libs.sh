#!/bin/bash

# Copy shared libraries
# (used internally - don't run this yourself!)

infiles=`find . -name "*.so"`
for infile in $infiles; do
  outfile=${DESTDIR}/usr/local/lib/pygeode/plugins/rpn/$infile
  mkdir -p ${outfile%/*}
  cp -i $infile $outfile
done
