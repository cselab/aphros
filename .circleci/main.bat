call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"
nmake /DCXXFLAGS="/O2 /nologo" /DCFLAGS="/O2 /nologo" /k /f Makefile_nmake
main.exe --logo --exit


