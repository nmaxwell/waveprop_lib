

icpc -Wall -fPIC -c waveprop.cpp -o libwaveprop.o $std_link
icpc -shared -Wl,-soname,libwaveprop.so.1 -o libwaveprop.so.1.0  libwaveprop.o  -lpng -lpngwriter -lz -lfreetype -lrt -pthread -lgsl -lgslcblas -larprec  -lfftw3_threads -lfftw3 -lz

#-lguide -pthread


cp libwaveprop.so.1.0 /usr/local/lib/
cp libwaveprop.so.1.0 /usr/lib/
cp libwaveprop.so.1.0 /lib64/
rm libwaveprop.o
rm libwaveprop.so.1.0
ln -sf /usr/local/lib/libwaveprop.so.1.0 /usr/local/lib/libwaveprop.so
ln -sf /usr/local/lib/libwaveprop.so.1.0 /usr/local/lib/libwaveprop.so.1
ln -sf /usr/lib/libwaveprop.so.1.0 /usr/lib/libwaveprop.so
ln -sf /usr/lib/libwaveprop.so.1.0 /usr/lib/libwaveprop.so.1
ln -sf /lib64/libwaveprop.so.1.0 /lib64/libwaveprop.so
ln -sf /lib64/libwaveprop.so.1.0 /lib64/libwaveprop.so.1