

icpc -Wall -fPIC -c waveprop.cpp -o libwaveprop.o  `freetype-config --cflags` -lpng -lpngwriter -lz -lfreetype -lrt -pthread -lgsl -lgslcblas -larprec  -lfftw3_threads -lfftw3 -I/usr/local/include -L/usr/local/lib -lz
icpc -shared -Wl,-soname,libwaveprop.so.1 -o libwaveprop.so.1.0  libwaveprop.o  `freetype-config --cflags` -lpng -lpngwriter -lz -lfreetype -lrt -pthread -lgsl -lgslcblas -larprec  -lfftw3_threads -lfftw3 -I/usr/local/include -L/usr/local/lib -lz


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


ln -svf  /workspace/waveprop_lib/python_waveprop/waveprop.py  /usr/lib/python2.6/site-packages/waveprop.py
