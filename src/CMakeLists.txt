include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/fftw)

set(SITUS_CPP matchpt.cpp voltrac.cpp volfltr.cpp)
set(AUX lib_err.c lib_eul.c lib_pio.c lib_pow.c lib_pwk.c lib_rnd.c lib_smp.c lib_std.c lib_tim.c lib_vec.c lib_vio.c lib_vwk.c lib_mpt.cpp lib_svt.cpp lib_sba.cpp)

add_library(situs SHARED ${AUX})
target_link_libraries(situs libfftw)

add_executable(map2map map2map.c)
target_link_libraries(map2map situs)

add_executable(volhist volhist.c)
target_link_libraries(volhist situs)

add_executable(voldiff voldiff.c)
target_link_libraries(voldiff situs)

add_executable(qvol qvol.c)
target_link_libraries(qvol situs)

add_executable(qpdb qpdb.c)
target_link_libraries(qpdb situs)

add_executable(pdb2sax pdb2sax.c)
target_link_libraries(pdb2sax situs)

add_executable(vol2pdb vol2pdb.c)
target_link_libraries(vol2pdb situs)

add_executable(colores colores.c)
target_link_libraries(colores situs)

add_executable(voledit voledit.c)
target_link_libraries(voledit situs)

add_executable(collage collage.c)
target_link_libraries(collage situs)

add_executable(pdbsymm pdbsymm.c)
target_link_libraries(pdbsymm situs)

add_executable(eul2pdb eul2pdb.c)
target_link_libraries(eul2pdb situs)

add_executable(qplasty qplasty.c)
target_link_libraries(qplasty situs)

add_executable(volaver volaver.c)
target_link_libraries(volaver situs)

add_executable(volmult volmult.c)
target_link_libraries(volmult situs)

add_executable(matchpt matchpt.cpp)
target_link_libraries(matchpt situs)

add_executable(voltrac voltrac.cpp)
target_link_libraries(voltrac situs)

add_executable(volfltr volfltr.cpp)
target_link_libraries(volfltr situs)
