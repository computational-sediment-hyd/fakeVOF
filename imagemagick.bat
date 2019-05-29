@echo off

set PATH=C:\Program Files\ImageMagick-7.0.8-Q16;%PATH%

magick convert -layers optimize -loop 0 -delay 5 *.png anim.gif
