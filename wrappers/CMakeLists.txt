enable_language( Fortran )
include( FortranCInterface )
FortranCInterface_HEADER( ${CMAKE_BINARY_DIR}/fortran_mangle.h
  MACRO_NAMESPACE "FC_"
  SYMBOL_NAMESPACE "FC_"
  )

list(APPEND includePath ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

find_package( BLAS   REQUIRED )
find_package( LAPACK   REQUIRED )
find_package( SCALAPACK REQUIRED)

add_library( blas_lapack_comm INTERFACE )

target_link_libraries( blas_lapack_comm INTERFACE
  BLAS::BLAS
  LAPACK::LAPACK
  SCALAPACK::SCALAPACK
)

target_include_directories( blas_lapack_comm INTERFACE
  "$<BUILD_INTERFACE:${includePath}>"
)

