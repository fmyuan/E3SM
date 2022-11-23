string(APPEND CPPDEFS " -DCPRWSL")
if (COMP_CLASS STREQUAL cpl)
  string(APPEND LDFLAGS " ")
endif()
