cmake_minimum_required(VERSION 2.8)

project(papara_nt)
set( CMAKE_NO_VERBOSE yes )

add_subdirectory(ublasJama-1.0.2.3)


# FIXME: the hack does not work, because IVYMIKE_INCLUDE can be cached...
# IF( NOT IVYMIKE_INCLUDE )
# SET( IVYMIKE_INCLUDE ./ivy_mike )

# this is a hack: only add subdirectories if IVYMIKE_INCLUDE is not already set.
# IVYMIKE_INCLUDE is abused as a flag indicating that papara should be build in 'library-only' mode...

if( NOT PAPARA_NO_IVYMIKE )
  add_subdirectory(ivy_mike)
endif()
# add_subdirectory(genassign_blast)
# ENDIF()

include_directories( ${include_directories} ${IVYMIKE_INCLUDE} . ublasJama-1.0.2.3 )



# handle the stupid autolinking shit of visualc++
IF(WIN32)

  IF( NOT BOOST_ROOT)
    set(BOOST_ROOT C:\\src\\boost_1_48_0)
  ENDIF()

  include_directories( ${BOOST_ROOT} )
  LINK_DIRECTORIES(C:/src/poco-1.4.5/lib )
  include_directories(C:/src/poco-1.4.5/Foundation/include) 

  set( BOOST_LIBS )
  set( SYSDEP_LIBS )
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_SCL_SECURE_NO_WARNINGS /INCREMENTAL:NO") # as long as there is no support for std::array, these warnings are plain stupid!

  file(GLOB ALL_HEADERS *.h)

ELSE(WIN32)

  include_directories( ${BOOST_ROOT} )
  if( USE_CPP11 )
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -march=native")
  else()
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++98 -pedantic -Wall -march=native")
  endif()
  set( BOOST_LIBS boost_thread boost_program_options)
  set(SYSDEP_LIBS pthread)
  #LINK_DIRECTORIES( ${LINK_DIRECTORIES} /usr/lib64/atlas-sse2 )

  set( ALL_HEADERS )
ENDIF(WIN32)



ADD_LIBRARY( papara_core STATIC papara.cpp pvec.cpp pars_align_seq.cpp pars_align_gapp_seq.cpp parsimony.cpp sequence_model.cpp align_utils.cpp blast_partassign.cpp )
set_property(TARGET papara_core PROPERTY CXX_STANDARD 11)

# add_executable(papara_nt main.cpp pvec.cpp pars_align_seq.cpp pars_align_gapp_seq.cpp parsimony.cpp ${ALL_HEADERS})
add_executable(papara papara2_main.cpp  ${ALL_HEADERS})
add_executable(fasta_random_sample2 fasta_random_sample2.cpp  ${ALL_HEADERS})
add_executable(fasta_to_phy fasta_to_phy.cpp  ${ALL_HEADERS})
add_executable(phy_to_fasta phy_to_fasta.cpp  ${ALL_HEADERS})


target_link_libraries(papara papara_core ivymike ublas_jama ${SYSDEP_LIBS} )
set_property(TARGET papara PROPERTY CXX_STANDARD 11)

target_link_libraries(phy_to_fasta ivymike ${SYSDEP_LIBS} )
set_property(TARGET phy_to_fasta PROPERTY CXX_STANDARD 11)

if(NOT WIN32)
# these two seem to imply a lib dependency to boost_date_time. NOPE, won't build that...
	add_executable(phy_cut_partition phy_cut_partition.cpp ${ALL_HEADERS})
	target_link_libraries(phy_cut_partition papara_core ublas_jama ivymike ${SYSDEP_LIBS} )
	set_property(TARGET phy_cut_partition PROPERTY CXX_STANDARD 11)

	add_executable(phy_megamerge phy_megamerge.cpp  ${ALL_HEADERS})
	target_link_libraries(phy_megamerge ivymike ${SYSDEP_LIBS} )
	set_property(TARGET phy_megamerge PROPERTY CXX_STANDARD 11)
endif()

#add_executable(stepwise_addition_pro stepwise_addition_pro.cpp pvec.cpp pairwise_seq_distance.cpp tree_similarity.cpp parsimony.cpp raxml_interface.cpp sequence_model.cpp ${ALL_HEADERS} )
#add_executable(propara propara.cpp pvec.cpp raxml_interface.cpp tree_similarity.cpp parsimony.cpp )
#add_executable(dump_anc_probs dump_anc_probs.cpp raxml_interface.cpp tree_similarity.cpp )



#target_link_libraries(stepwise_addition_pro ivymike ${BOOST_LIBS} ${SYSDEP_LIBS} PocoFoundation ublas_jama)
#set_property(TARGET stepwise_addition_pro PROPERTY CXX_STANDARD 11)
#target_link_libraries(propara papara_core ivymike ${BOOST_LIBS} ${SYSDEP_LIBS} PocoFoundation ublas_jama)
#set_property(TARGET propara PROPERTY CXX_STANDARD 11)
#target_link_libraries(dump_anc_probs ivymike ublas_jama ${BOOST_LIBS} ${SYSDEP_LIBS} PocoFoundation)
#set_property(TARGET dump_anc_probs PROPERTY CXX_STANDARD 11)

install_targets( /bin papara )
