# Working with Wavelet Libraries #
## Wavelet Libraries in Win32 Environment ##
### Microsoft Visual C++ ###

This package contains two DLLs for MSVC++ - one Debug and one Release
version.The list of required files in order to compile wavelet2d based programs-

#Header File- wavelet2d.h

#Import Library- wavelet2d.lib

#Wavelet DLL - wavelet2d.dll

#FFTW DLL - libfftw3-3.dll

Header file is in the source (’src’) folder. Wavelet Libraries are in the respective Debug and Release folders. FFTW DLL is in the FFTW3 folder. Alternatively, you may chose to install FFT library from www.fftw.org. The FFTW sourcecodes are also available at FFTW website under GNU-GPL license.

In order to use wavelet libraries, the easiest way is to add import library
(wavelet2d.lib) to additional dependency in your project which in turn will
handle the DLLs for you. The two DLLs must be in your program path.If you
are new to MSVC, you may want to learn more about DLL search path at

[DLL Search Path](http://msdn.microsoft.com/en-us/library/7d83bc18(VS.71).aspx)

You may also chose to directly call DLLs from your program.Regardless
of how you handle DLLs, mixing debug and release versions is not a good
idea.Both versions are named wavelet2d.dll so putting them in system path is
not a good idea either. Preferred way is to use Visual Studio project settings
to link to the wavelet2d.lib file and put the release and debug wavelet2d dlls
along with libfftw-3.3dll in the directory containing the respective release and
debug executables of your project. Import library will open the DLLs.

### GNU-GCC (MinGW) Based Compilation ###

If you use GCC compiler either through MSYS or from one of the IDEs(Codeblocks,
Eclipse, QTcreator, Netbeans etc.),I have also added both debug and release
versions of wavelet2d DLLs and wavelet2s static libraries. Required Files for
MinGW DLLS-

• Header Files- wavelet2d.h

• Import Library- libwavelet2d.dll.a

• Wavelet DLL - wavelet2d.dll

• FFTW DLL - libfftw3-3.dll (Import library is also included)

It is recommended that you link to the Release and Debug folders as they are.
Once you link to import library and specify the folder location, wavelet import
library should automatically open the two DLLs. Mixing Debug and Release
versions shouldn’t pose any problems in this case.

Working with MinGW Static libraries is even more straightforward. You
need to include wavelet2s.h header in your program , link to static .a library and you are set.

**Known Issue - Some newer versions of MinGW-g++ compiler have problem
working with .o codes(and hence static libraries) compiled in older versions of MinGW (3.x). If you are facing this problem try MINGW4.5 compiled static libraries or build wavelet2d from the source code with your own compiler. You'll need to statically build fftw3 library first if you are going the source code way. The DLLs should be working fine in any case so that's always another option.**


## Working in LINUX Environment ##

1. This package contains two wavelet libraries- libwavelet2d.so.1.0 (shared)
and libwavelet2s.a (static) compiled essentially from the same source code.
Source code is available in the ’src’ folder.

2. You may need to link to header files that are included with their respective libraries. They are also available in the ’src’ folder.

3. You may want to install shared library in one of your existing paths to
make compilation easier. You can create sym links once inside the folder(
eg., /usr/local/lib where you are installing libwavelet2d.so.1.0) by using following commands

> ln -sf libwavelet2d.so.1.0 libwavelet2d.so

> ln -sf libwavelet2d.so.1.0 libwavelet2d.so.1

> ldconfig

You will probably need superuser privileges to perform previous steps. If you
don’t have su privilege then you can move libwavelet2d.so.1.0 to your work
folder or where your source code is, create sym links as before and then put
your folder in the path during runtime.

ln -sf libwavelet2d.so.1.0 libwavelet2d.so

ln -sf libwavelet2d.so.1.0 libwavelet2d.so.1

export LD\_LIBRARY\_PATH=.


libwavelet2s.a :: Working with static library is pretty straightforward. You will
only need to include wavelet2s.h in your program and specify the library
(-lwavelet2s flag) , include path (-I<path to wavelet2s.h> flag) and library
path (-L<path to libwavelet2s.a> flag).