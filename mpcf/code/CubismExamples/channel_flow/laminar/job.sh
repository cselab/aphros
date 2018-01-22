#!/bin/bash
PATH="$(pwd -P)":$PATH
export OMP_NUM_THREADS=12

mpcf-node -simConf channelFlow.conf "$@"
code=$?
# mail -s "$(hostname): Job \"channel_node\" has finished" fabianw@mavt.ethz.ch < EOT

exit $code
