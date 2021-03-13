#!/bin/sh -eu

cd $REPO/deploy
./install_setenv $HOME/.local/aphros
. $HOME/.local/bin/ap.setenv
make -j8
make install

cd $REPO/src
make -j8
