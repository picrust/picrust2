set( NACL ON )

set( PLATFORM_TRIPLET           "i686-nacl" )
set( PLATFORM_PREFIX            "/space/src/nacl_sdk/pepper_23/toolchain/linux_x86_newlib" )

set( CMAKE_SYSTEM_NAME          "Linux" )
set( CMAKE_SYSTEM_PROCESSOR     "i686" )
set( CMAKE_FIND_ROOT_PATH       "${PLATFORM_PREFIX}/${PLATFORM_TRIPLET}" )
set( CMAKE_C_COMPILER           "${PLATFORM_PREFIX}/bin/${PLATFORM_TRIPLET}-gcc" )
set( CMAKE_CXX_COMPILER         "${PLATFORM_PREFIX}/bin/${PLATFORM_TRIPLET}-g++" )
set( CMAKE_STRIP                "${PLATFORM_PREFIX}/bin/${PLATFORM_TRIPLET}-strip" )

set( LUA_MATH_LIBRARY           "${PLATFORM_PREFIX}/x86_64-nacl/lib32/libm.a" CACHE PATH "" )
set( PEPPER_LIBRARY             "${PLATFORM_PREFIX}/x86_64-nacl/lib32/libppapi.a" CACHE PATH "" )
set( PEPPER_CXX_LIBRARY         "${PLATFORM_PREFIX}/x86_64-nacl/lib32/libppapi_cpp.a" CACHE PATH "" )
set( PEPPER_GLES2_LIBRARY       "${PLATFORM_PREFIX}/x86_64-nacl/lib32/libppapi_gles2.a" CACHE PATH "" )

