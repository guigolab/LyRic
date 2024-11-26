#!/bin/bash


#while true; do
#date;
echo "Forced unmount just in case..."
fusermount -uz /nfs
echo "Mounting..."
#sshfs -o reconnect -C -o workaround=all  jlagarde@ant-login.linux.crg.es:/nfs/ /nfs/
sshfs -o reconnect,compression=yes,auto_cache,cache_timeout=5,transform_symlinks,ServerAliveInterval=60,ServerAliveCountMax=3,ssh_command='autossh -M 0' jlagarde@ant-login.linux.crg.es:/nfs/ /nfs/

#su - autossh -s /bin/bash -c "/usr/bin/sshfs -o reconnect,compression=yes,auto_cache,cache_timeout=5,transform_symlinks,allow_other,idmap=user,ServerAliveInterval=60,ServerAliveCountMax=3,ssh_command='autossh -M 0' autossh@MASTER_SERVER_IP:/ /opt/sshfs -p MASTER_SERVER_SSH_PORT"  #echo "Done"


#    sleep 5m
#done &> ~jlagarde/mountsshfsCrg.sh.log
