#!/bin/bash

set -o verbose

now=`date +"%m_%d_%Y"`
email="zacharyanderson2024@u.northwestern.edu"
backlog="/home/zaz3744/backup_log_$now.txt"

#rm .netrc

yes "" | openssl enc -in .netrc.enc \
    -d -aes-256-cbc \
    -pass stdin > .netrc

lftp ftp.box.com <<Backup

lcd /home/zaz3744/ACNlab/projects/BrainMAPD_func_conn
cd Quest_Backup
mirror --reverse --verbose --only-newer --log="$backlog" #-c essentially does rsync
quit

# Backup

cd ~
# rm .netrc

mail -s "backup_log_$now" "$email" < "$backlog"
