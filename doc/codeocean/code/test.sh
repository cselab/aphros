#!/bin/sh

set -eu

REPO=/code/aphros

cd $REPO/deploy
./install_setenv $HOME/.local/aphros
. $HOME/.local/bin/ap.setenv
make -j8
make install

cd $REPO/src
make -j8

cd $REPO/src
make test || true

rsync -av $REPO/src/build/Testing /results/