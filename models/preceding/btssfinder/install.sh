INSTALLATION

#1. check for gfortran
#if not installed, wget tarball (.tar.gz) and install gfortran
#2. Add gfortran to USER $PATH - e.g. ~/.bash_profile
#3. Check if fortran works
#3. Identify the current location of BTSFINDER and Add bTSFINDER_data to user $PATH - e.g. ~/.bash_profile

CONDITION=which gfortran
if $CONDITION ==0; then
	echo "gfotran was not detected on this machine. Downloading the file from echo http://gfortran.meteodat.ch/download/x86_64/gcc-4.8-infrastructure.tar.xz"
	wget http://gfortran.meteodat.ch/download/x86_64/gcc-4.8-infrastructure.tar.xz
fi
tar -xvf gcc-4.8-infrastructure.tar.xz

sed -i -r 'PATH=~/fortran/:$PATH' ~/.bash_profile
source ~/.bash_profile
fortran --version
wget ''
