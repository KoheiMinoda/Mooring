# ! /bin/bash

cd mbdyn
rm -f *.o
cd ../

cd modules/module-aerodyn13_b1.09
rm -rf .libs
rm -f *.la *.lo *.o
cd source
chmod u+x ./compileAD.sh
./compileAD.sh
cd ../../../

# build modules
LDFLAGS=-rdynamic
LIBS=/usr/lib/x86_64-linux-gnu/libltdl.a
CC="gcc-9" CXX="g++-9" F77="gfortran-9" FC="gfortran-9" CPPFLAGS=-I/usr/include/suitesparse ./configure --enable-runtime-loading 
CC="gcc-9" CXX="g++-9" F77="gfortran-9" FC="gfortran-9" CPPFLAGS=-I/usr/include/suitesparse ./configure --with-module="aerodyn13_b1.09"
make
sudo make install

# delete .mod file for viewing
cd modules/module-aerodyn13_b1.09
rm -f *.mod
cd ../

# install ufmpack
CPPFLAGS=-I/usr/include/suitesparse ./configure
make
sudo make install
