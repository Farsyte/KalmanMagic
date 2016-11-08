#!/bin/sh -
set -e

echo
echo 'Cyclomatic Complexity:'
pmccabe -v /dev/null
pmccabe "$@" | grep -v ': mexFunction$$'
echo
echo 'Source Line Counts:'
sloccount --details "$@" > /dev/null 2>&1
sloccount --cached --details "$@"
echo
echo 'Baseline COCOMO (*unconfigured*):'
sloccount --cached           "$@"
echo
