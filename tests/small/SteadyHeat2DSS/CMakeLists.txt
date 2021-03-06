# Name the test with the directory name
get_filename_component(testNameRoot ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# Copy mesh files from source to binary dir
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quadQuadSS.exo
               ${CMAKE_CURRENT_BINARY_DIR}/quadQuadSS.exo COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quadQuadSS.exo.4.0
               ${CMAKE_CURRENT_BINARY_DIR}/quadQuadSS.exo.4.0 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quadQuadSS.exo.4.1
               ${CMAKE_CURRENT_BINARY_DIR}/quadQuadSS.exo.4.1 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quadQuadSS.exo.4.2
               ${CMAKE_CURRENT_BINARY_DIR}/quadQuadSS.exo.4.2 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quadQuadSS.exo.4.3
               ${CMAKE_CURRENT_BINARY_DIR}/quadQuadSS.exo.4.3 COPYONLY)

####################################
###          Main tests          ###
####################################

set(testName ${testNameRoot})

if (ALBANY_EPETRA) 
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)

  add_test(${testName}_Epetra ${Albany.exe} input.yaml)
  set_tests_properties(${testName}_Epetra PROPERTIES LABELS "Basic;Epetra;Forward")
endif()

if (ALBANY_IFPACK2)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/inputT.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/inputT.yaml COPYONLY)

  add_test(${testName}_Tpetra ${Albany.exe} inputT.yaml)
  set_tests_properties(${testName}_Tpetra PROPERTIES LABELS "Basic;Tpetra;Forward")
endif ()

####################################
###          DuDx tests          ###
####################################

set(testName ${testNameRoot}_dydx)

if (ALBANY_EPETRA) 
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_dudx.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_dudx.yaml COPYONLY)

  add_test(${testName}_Epetra ${Albany.exe} input_dudx.yaml)
  set_tests_properties(${testName}_Epetra PROPERTIES LABELS "Basic;Epetra;Forward")
endif()

if (ALBANY_IFPACK2)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/inputT_dudx.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/inputT_dudx.yaml COPYONLY)

  add_test(${testName}_Tpetra ${Albany.exe} inputT_dudx.yaml)
  set_tests_properties(${testName}_Tpetra PROPERTIES LABELS "Basic;Tpetra;Forward")
endif ()
