#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# Take an optional argument: the basename of the game file.
game="maze2"
if [ "$#" -ge 1 ]; then
  game=$1
fi

Player1=UnClonRandom
Player2=Dummy	
Player3=UnDummyMes
Player4=Sanfe

VIEWER_PATH=Viewer

./Game $Player1 $Player2 $Player3 $Player4 -s 34534 < ${game}.cnf > $VIEWER_PATH/${game}.res

start firefox "file:///C:/Users/Gonzalo/Documents/Programacion/C++_Workspace_Qt-Creator/StarWar/$VIEWER_PATH/viewer.html?game=${game}.res"
