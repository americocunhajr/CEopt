! /bin/bash
 
## look in comments
## http://opendevice.blogspot.com/2007/05/eps-to-pdf-how-to-avoid-clipping.html
for file in "$@"
do
ps2pdf -dEPSCrop "$file"
done