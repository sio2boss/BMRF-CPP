#! /bin/sh

if [ -e /vagrant/devops/provision.stop ]; then
   echo "Shell Provision Bypassed"
   exit
fi


echo "Shell Provision: Update aptitude package cache"
sleep 30
sudo aptitude update

echo "Shell Provision: Install developer utils"
sudo aptitude -y install vim ctags subversion git-core curl screen gdb patch build-essential node zlibc cmake python unzip liblog4cxx10-dev libboost-all-dev 

echo "Shell Provision: Install threadpool"
cd /tmp
wget -q http://downloads.sourceforge.net/project/threadpool/threadpool/0.2.5%20%28Stable%29/threadpool-0_2_5-src.zip
unzip threadpool-0_2_5-src.zip  &> /dev/null
sudo cp -r threadpool-0_2_5-src/threadpool/boost/* /usr/include/boost/
rm -rf /tmp/threadpool-0_2_5-src*

echo "Shell Provision: Install HDF5"
cd /tmp
wget -q http://www.hdfgroup.org/ftp/HDF5/current/bin/linux-x86_64/hdf5-1.8.11-linux-x86_64-shared.tar.gz &> /dev/null
tar xvfz hdf5-1.8.11-linux-x86_64-shared.tar.gz  &> /dev/null
sudo cp -r hdf5-1.8.11-linux-x86_64-shared/include/* /usr/local/include/
sudo cp -r hdf5-1.8.11-linux-x86_64-shared/bin/*     /usr/local/bin/
sudo cp -r hdf5-1.8.11-linux-x86_64-shared/lib/*     /usr/local/lib/
sudo cp -r hdf5-1.8.11-linux-x86_64-shared/share/*   /usr/local/share/
rm -rf hdf5-1.8.11-linux-x86_64-shared*
cd /usr/local/bin/
sudo ./h5redeploy -prefix=/usr/local -force

### Add stuff here ###

echo "Shell Provision: Updating locate cache"
sudo updatedb

echo "Shell Provision: Done"
touch /vagrant/devops/provision.stop
