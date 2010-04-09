
g++ -Wall -fPIC -c waveprop.cc -o libwaveprop.o $std_link
g++ -shared -Wl,-soname,libwaveprop.so.1 -o libwaveprop.so.1.0   libwaveprop.o $std_link
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


