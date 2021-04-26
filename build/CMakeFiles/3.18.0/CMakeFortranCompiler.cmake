set(CMAKE_Fortran_COMPILER "/p/software/juwelsbooster/stages/2020/software/GCCcore/9.3.0/bin/gfortran")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "9.3.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/p/software/juwelsbooster/stages/2020/software/binutils/2.34-GCCcore-9.3.0/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/p/software/juwelsbooster/stages/2020/software/GCCcore/9.3.0/bin/gcc-ar")
set(CMAKE_RANLIB "/p/software/juwelsbooster/stages/2020/software/binutils/2.34-GCCcore-9.3.0/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "/p/software/juwelsbooster/stages/2020/software/GCCcore/9.3.0/bin/gcc-ranlib")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/p/software/juwelsbooster/stages/2020/software/GCCcore/9.3.0/lib/gcc/x86_64-pc-linux-gnu/9.3.0/finclude;/p/software/juwelsbooster/stages/2020/software/Boost/1.74.0-gpsmpi-2020/include;/p/software/juwelsbooster/stages/2020/software/ICU/67.1-GCCcore-9.3.0/include;/p/software/juwelsbooster/stages/2020/software/bzip2/1.0.8-GCCcore-9.3.0/include;/p/software/juwelsbooster/stages/2020/software/ncurses/6.2-GCCcore-9.3.0/include;/p/software/juwelsbooster/stages/2020/software/imkl/2020.2.254-gpsmpi-2020/mkl/include/fftw;/p/software/juwelsbooster/stages/2020/software/imkl/2020.2.254-gpsmpi-2020/mkl/include;/p/software/juwelsbooster/stages/2020/software/psmpi/5.4.7-1-GCC-9.3.0/include;/p/software/juwelsbooster/stages/2020/software/libxml2/2.9.10-GCCcore-9.3.0/include/libxml2;/p/software/juwelsbooster/stages/2020/software/libxml2/2.9.10-GCCcore-9.3.0/include;/p/software/juwelsbooster/stages/2020/software/XZ/5.2.5-GCCcore-9.3.0/include;/p/software/juwelsbooster/stages/2020/software/pscom/5.4-default/include;/p/software/juwelsbooster/stages/2020/software/UCX/1.9.0/include;/p/software/juwelsbooster/stages/2020/software/CUDA/11.0/nvvm/include;/p/software/juwelsbooster/stages/2020/software/CUDA/11.0/extras/CUPTI/include;/p/software/juwelsbooster/stages/2020/software/CUDA/11.0/include;/p/software/juwelsbooster/stages/2020/software/numactl/2.0.13/include;/p/software/juwelsbooster/stages/2020/software/binutils/2.34-GCCcore-9.3.0/include;/p/software/juwelsbooster/stages/2020/software/zlib/1.2.11-GCCcore-9.3.0/include;/p/software/juwelsbooster/stages/2020/software/GCCcore/9.3.0/lib/gcc/x86_64-pc-linux-gnu/9.3.0/include;/p/software/juwelsbooster/stages/2020/software/GCCcore/9.3.0/include;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;gcc_s;gcc;quadmath;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/p/software/juwelsbooster/stages/2020/software/psmpi/5.4.7-1-GCC-9.3.0/lib64;/p/software/juwelsbooster/stages/2020/software/pscom/5.4-default/lib64;/p/software/juwelsbooster/stages/2020/software/UCX/1.9.0/lib64;/p/software/juwelsbooster/stages/2020/software/nvidia-driver/default/lib64;/p/software/juwelsbooster/stages/2020/software/GCCcore/9.3.0/lib/gcc/x86_64-pc-linux-gnu/9.3.0;/p/software/juwelsbooster/stages/2020/software/GCCcore/9.3.0/lib64;/lib64;/usr/lib64;/p/software/juwelsbooster/stages/2020/software/Boost/1.74.0-gpsmpi-2020/lib;/p/software/juwelsbooster/stages/2020/software/ICU/67.1-GCCcore-9.3.0/lib;/p/software/juwelsbooster/stages/2020/software/bzip2/1.0.8-GCCcore-9.3.0/lib;/p/software/juwelsbooster/stages/2020/software/ncurses/6.2-GCCcore-9.3.0/lib;/p/software/juwelsbooster/stages/2020/software/imkl/2020.2.254-gpsmpi-2020/mkl/lib/intel64;/p/software/juwelsbooster/stages/2020/software/imkl/2020.2.254-gpsmpi-2020/lib/intel64;/p/software/juwelsbooster/stages/2020/software/psmpi/5.4.7-1-GCC-9.3.0/lib;/p/software/juwelsbooster/stages/2020/software/libxml2/2.9.10-GCCcore-9.3.0/lib;/p/software/juwelsbooster/stages/2020/software/XZ/5.2.5-GCCcore-9.3.0/lib;/p/software/juwelsbooster/stages/2020/software/UCX/1.9.0/lib;/p/software/juwelsbooster/stages/2020/software/CUDA/11.0/lib64/stubs;/p/software/juwelsbooster/stages/2020/software/CUDA/11.0/lib64;/p/software/juwelsbooster/stages/2020/software/numactl/2.0.13/lib;/p/software/juwelsbooster/stages/2020/software/binutils/2.34-GCCcore-9.3.0/lib;/p/software/juwelsbooster/stages/2020/software/zlib/1.2.11-GCCcore-9.3.0/lib;/p/software/juwelsbooster/stages/2020/software/GCCcore/9.3.0/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
