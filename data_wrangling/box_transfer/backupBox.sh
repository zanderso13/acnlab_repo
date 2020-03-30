#!/bin/bash

set -o verbose

now=`date +"%m_%d_%Y"`
email="btklm07@gmail.com"
backlog="/home/btk2142/backup_log_$now.txt"

#rm .netrc

yes "" | openssl enc -in .netrc.enc \
    -d -aes-256-cbc \
    -pass stdin > .netrc

lftp ftp.box.com <<Backup

lcd /projects/b1081/
cd Quest_Backup
mirror --reverse --verbose --only-newer --log="$backlog"
quit

Backup

cd ~
rm .netrc

mail -s "backup_log_$now" "$email" < "$backlog"
