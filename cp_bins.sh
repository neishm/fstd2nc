#!/bin/bash

# Copy binaries to installation destination
# (used internally - don't run this yourself!)

bindir=${DESTDIR}/usr/local/bin
mkdir -p $bindir
cp ./bin/* ${bindir}/

