#!/bin/sh

rm pngs/frame00000.png

mencoder mf://pngs/*.png -mf fps=24 -ovc lavc -lavcopts vcodec=msmpeg4v2 -oac copy -o output.avi

ffmpeg -i "output.avi" -acodec libmp3lame -ab 192 "output.mov"

rm output.avi
