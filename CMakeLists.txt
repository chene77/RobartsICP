CMAKE_MINIMUM_REQUIRED( VERSION 2.2 )

PROJECT( RobartsICP )


INCLUDE_DIRECTORIES(
  C:/chene/src/ann_1.1.2/include
  )

LINK_DIRECTORIES(
  C:/chene/src/ann_1.1.2/MS_Win32/dll/Release
  )

ADD_EXECUTABLE( ASICP
  asicp_main.cpp
  mathUtils.cpp
  mathUtils.h
  matrix.h
  jacobi.h
  svd.h
  ASOPP_Major.cpp
  ASOPP_Major.h
  pointRegistration.cpp
  pointRegistration.h
  asicp.cpp
  asicp.h
)
TARGET_LINK_LIBRARIES( ASICP
  ANN
)


ADD_EXECUTABLE( icp
  icp_main.cpp
  icp.cpp
  icp.h
  pointRegistration.cpp
  pointRegistration.h
  mathUtils.cpp
  mathUtils.h
  svd.h
  matrix.h
)

TARGET_LINK_LIBRARIES( icp
  ANN
)

ADD_EXECUTABLE( robust_icp
  robust_icp_main.cpp
  robust_icp.cpp
  robust_icp.h
  tukey_pointRegistration.h
  tukey_pointRegistration.cpp
  heterogeneous_pointRegistration.h
  heterogeneous_pointRegistration.cpp
  mathUtils.h
  mathUtils.cpp
  qsort.h
)
TARGET_LINK_LIBRARIES( robust_icp
  ANN
  )
