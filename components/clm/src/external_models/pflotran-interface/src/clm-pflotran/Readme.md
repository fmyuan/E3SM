# INSTRUCTIONS #

This is the clm-pflotran interface, including minimal PFLOTRAN codes as an external model, for coupling PFLOTRAN into CLM in DOE sponsored NGEE-Arctic Project. 

The model coupling aims to provide a full alternative solution for CLM-CN's subsurface C/N biogeochemistry and thermal-hydrology, i.e. PFLOTRAN.

***This interface is specifically released for ```PETSC v.3.9x or above```, and ```E3SM v1.2 or above```. ***


## How do I get set up? ##

**(1)** *git clone the repository*, if not yet properly updating codes. (OPTIONAL)
```
git clone https://github.com/fmyuan/E3SM.git
```

*NOTE: when git clone, checkout branch 'elm-pflotran-II'*


**(2)** *if coupling CLM with PFLOTRAN, need to soft-link required original PFLOTRAN source codes.* **NOTE: ELM will be compiled including full PFLOTRAN interface source codes AUTOMATICALLY**
```
cd ./compoenents/clm/src/external_models/pflotran-interface/src/clm-pflotran
make link_common_src

(NOTE: if really don't want your codes large, you could reverse the link, by issue command:
make clean_common_src)
```

***(3) FINALLY***, build ELM as usual, BUT must do modifying ELM's configuration .xml as following.

*I.* **PETSc library:** You must be sure that $PETSC_PATH are defined in your bash environment setting. AND edit **'env_build.xml'** as following:
```
  <group id="build_component_clm">
    <entry id="CLM_USE_PETSC" value="TRUE">
```

*II.* **env_mach_specific.xml** editing for each supported machine (example CADES at ORNL) after './case.setup', to turn on relevant options.

```
   <environment_variables>
     <env name="CLM_PFLOTRAN_COUPLED">TRUE</env>
     <env name="CLM_PFLOTRAN_COLMODE">TRUE</env>
     <!-- by blanking the following 2 names, PETSC libs excluded in e3sm.exe when NOT coupling with PFLOTRAN -->
     <env name="PFLOTRAN_INC"> -I$ENV{PETSC_DIR}/include</env>
     <env name="PFLOTRAN_LIB"> -L$ENV{PETSC_DIR}/lib -lpetsc -lmetis -lparmetis</env>
   </environment_variables>

```
**NOTE:** You DON'T need to set $CLM_PFLOTRAN_SOURCE_DIR in your bash environment setting. If you DON'T want to include PETSc library in your e3sm.exe (and of course no use of coupled models ), edit 'env_mach_specific.xml' as following:
```
<environment_variables>
  <env name="CLM_PFLOTRAN_COUPLED">FALSE</env>
  <env name="CLM_PFLOTRAN_COLMODE">FALSE</env>
<!--
  <env name="PFLOTRAN_INC"> -I$ENV{PETSC_PATH}/include</env>
  <env name="PFLOTRAN_LIB"> -L$ENV{PETSC_PATH}/lib -lpetsc -lmetis -lparmetis</env>
-->
</environment_variables>
```

*III.* **env_run.xml** editing to turn on RUN TIME option. **NOTE:** you can turn it off, with coupled codes fully compiled, so that same model for either coupled run or not.

```
    <entry id="CLM_INTERFACE_MODE" value="pflotran">
      <type>char</type>
      <valid_values>off,bgc,pflotran</valid_values>
      <desc>CLM build-namelist options for ' -clm_interface_mode'
         FOR run-time options to apply for new CLM capabilities via clm-interface.
         VALID options: (1) off (default);
                        (2) bgc - run native clm BGC via interface (for purpose of testing interface);
                        (3) pflotran - run coupled clm-pflotran via interface (coupled module(s) upon pflotran.in).</desc>
    </entry>

```


***UPDATED: 2019-05-16***
