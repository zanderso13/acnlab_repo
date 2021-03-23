#!/bin/bash

cd ~

openssl enc -in .netrc.enc \
    -d -aes-256-cbc \
    -pass stdin > .netrc

set -e

inputfilepath="$1"
outputfilepath="$2"
#inputfolder="$3"
#filename="$4"

#ftp -i ftp.box.com

lftp ftp.box.com <<EOF

lcd "$outputfilepath"
cd "$inputfilepath"
mget *


quit

EOF

#mget "$filename"
#mget -c "$inputfolder"

cd ~
rm .netrc
