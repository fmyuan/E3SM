add_subdirectory(base)
add_subdirectory(betr)
add_subdirectory(fates)
add_subdirectory(ptm)
add_subdirectory(vsfm)
add_subdirectory(stub)
add_subdirectory(alquimia)
add_subdirectory(ats)

set(EMI_LIBRARIES emi_em_base; emi_em_betr; emi_em_fates; emi_em_ptm; emi_em_vsfm; emi_em_stub; emi_em_alquimia; emi_em_ats; ${EMI_LIBRARIES} PARENT_SCOPE)


