# here we piggy-back on EMI to get this installed
set(EM_ATS_SOURCES
  ExternalModelATS_readnlMod.F90
  ExternalModelATS_variables.F90
  ExternalModelATS_interface.F90
  ExternalModelATS.F90
)

# work toward replacing with compset and
# configure scripts for ELM. Possible
# that much of this can be deleted in the end.
include(add_emi_library)
add_emi_library(em_ats ${EM_ATS_SOURCES})

set(EMI_LIBRARIES em_ats;${EMI_LIBRARIES} PARENT_SCOPE)
set(EMI_LIBRARIES em_ats;${EMI_LIBRARIES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS em_ats DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()





