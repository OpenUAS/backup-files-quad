cmake_minimum_required(VERSION 2.8)

project(quad_ipopt)

SET(MANDATORY_LIBRARIES "")

#We look for ipopt libs and includes
FIND_LIBRARY(IPOPT_LIBRARY NAMES ipopt)
FIND_PATH(IPOPT_INCLUDE_DIR NAMES coin/IpTNLP.hpp)
MESSAGE(STATUS "IPOPT library: ${IPOPT_LIBRARY}")
MESSAGE(STATUS "IPOPT include dir: ${IPOPT_INCLUDE_DIR}")
IF(NOT IPOPT_LIBRARY OR NOT IPOPT_INCLUDE_DIR)
	MESSAGE(FATAL_ERROR "IPOPT support was requested, but ipopt library and headers could not be located. Please make sure that the linker can find the ipopt library and headers. The headers are assumed to be installed in a in coin/ directory")
ENDIF(NOT IPOPT_LIBRARY OR NOT IPOPT_INCLUDE_DIR)

#We look for lapack
FIND_PACKAGE(LAPACK REQUIRED)
MESSAGE(STATUS "lapack libraries: ${LAPACK_LIBRARIES}")
MESSAGE(STATUS "lapack linker flags: ${LAPACK_LINKER_FLAGS}")

#We look for dl
FIND_LIBRARY(DL_LIBRARY NAMES dl)
MESSAGE(STATUS "dl library: ${DL_LIBRARY}")
IF(NOT DL_LIBRARY)
	MESSAGE(FATAL_ERROR "IPOPT support was requested, but dl library could not be found.")
ENDIF(NOT DL_LIBRARY)

#We look for coinhsl
FIND_LIBRARY(COINHSL_LIBRARY NAMES coinhsl)
IF(COINHSL_LIBRARY)
	MESSAGE(STATUS "COINHSL Library Found ... linking it in")
	MESSAGE(STATUS "coin hsl library: ${COINHSL_LIBRARY}")
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${COINHSL_LIBRARY})
ENDIF(COINHSL_LIBRARY)

#We look for coinmetis
FIND_LIBRARY(COINMETIS_LIBRARY NAMES coinmetis)
IF(COINMETIS_LIBRARY)
	MESSAGE(STATUS "COINMETIS Library Found ... linking it in")
	MESSAGE(STATUS "coin metis library: ${COINMETIS_LIBRARY}")
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${COINMETIS_LIBRARY})
ENDIF(COINMETIS_LIBRARY)

#We look for coinmumps
FIND_LIBRARY(COINMUMPS_LIBRARY NAMES coinmumps)
IF(COINMUMPS_LIBRARY)
	MESSAGE(STATUS "COINMUMPS Library Found ... linking it in")
	MESSAGE(STATUS "coin mumps library: ${COINMUMPS_LIBRARY}")
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${COINMUMPS_LIBRARY})
ENDIF(COINMUMPS_LIBRARY)

#We check we have at least a linear solver
IF(NOT COINHSL_LIBRARY AND NOT COINMUMPS_LIBRARY)
	MESSAGE(STATUS "Linear Solver for IPOPT could not be detected. I tried with HSL and MUMPS. You need to link manually the library if necessary")
ENDIF(NOT COINHSL_LIBRARY AND NOT COINMUMPS_LIBRARY)
	
#We link all in
SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${IPOPT_LIBRARY} ${LAPACK_LIBRARIES} ${DL_LIBRARY} ${COINHSL_LIBRARY})
INCLUDE_DIRECTORIES("${IPOPT_INCLUDE_DIR}")

SET(QUAD_IPOPT_LIB_SRC_LIST
	${CMAKE_CURRENT_SOURCE_DIR}/quadNLP.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/quadFunctions.cpp
)
ADD_LIBRARY(quad_ipopt_static STATIC ${QUAD_IPOPT_LIB_SRC_LIST})
ADD_EXECUTABLE(quad_ipopt quadNLPmain.cpp)
TARGET_LINK_LIBRARIES(quad_ipopt ${MANDATORY_LIBRARIES} quad_ipopt_static)
