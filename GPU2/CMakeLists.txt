CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(nbody)

ENABLE_LANGUAGE(Fortran)
ENABLE_LANGUAGE(CXX)

set(CMAKE_CXX_STANDARD 11)
if(${CMAKE_VERSION} VERSION_GREATER "3.12.0")
cmake_policy(SET CMP0074 NEW)
endif()

option(USE_CUB
    "On to use CUB sorting"
  OFF
  )
set(USE_MPI OFF CACHE BOOL "Use of MPI is disabled." FORCE)

FIND_PACKAGE(CUDA REQUIRED)
FIND_PACKAGE(OpenMP REQUIRED)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} ${OpenMP_FORTRAN_FLAGS}")

OPTION(USE_ARCHAIN "Select the ARChain mode" ON)
OPTION(USE_KSBLOCK "Use block KS regularisation" OFF)
OPTION(SANITIZE "Inject address sanitization code" OFF)
OPTION(PROFILER "Link with profiler library" OFF)
OPTION(TRAP_NAN "Trap invalid floating point" OFF)

set (COMPILER_FLAGS)
if(SANITIZE)
  set (COMPILER_FLAGS "${COMPILER_FLAGS} -fsanitize=address")
endif()
if (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
  set (COMPILER_FLAGS "${COMPILER_FLAGS} -msse4 -march=native -mavx")
endif()

set (COMPILER_FLAGS
  "${COMPILER_FLAGS} -Wno-strict-aliasing")
set (COMPILER_FLAGS
  "${COMPILER_FLAGS} -Wno-unused-variable -Wno-unused-function")
set (COMPILER_FLAGS
  "${COMPILER_FLAGS} -fPIC")

set(MCMODEL "large" CACHE STRING "The mcmodel flag to use.")
if (MCMODEL)
  set (COMPILER_FLAGS
    "${COMPILER_FLAGS} -mcmodel=${MCMODEL}")
endif()
ADD_DEFINITIONS(-DHUGE=1.0e20 -DWITH_CUDA5 -D_FORCE_INLINES -DCUB_IGNORE_DEPRECATED_CPP_DIALECT)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reorder")

SET(DEBUG_FLAGS "-fstack-protector -fno-omit-frame-pointer -O0 -g")
set(CMAKE_EXE_LINKER_FLAGS
  "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")


set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${COMPILER_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILER_FLAGS}")
SET(CMAKE_CXX_FLAGS_DEBUG  "${DEBUG_FLAGS}")
SET(CMAKE_Fortran_FLAGS_DEBUG "-fcheck=all ${DEBUG_FLAGS}")
SET(CMAKE_CXX_FLAGS_RELEASE  "-O3")
SET(CMAKE_Fortran_FLAGS_RELEASE "-O3")
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--no-relax")

if (TRAP_NAN)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffpe-trap=invalid")
endif()

if (PROFILER)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ggdb")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lprofiler")
endif()

set(CUDA_OPTIONS
    -arch=sm_50
    -gencode=arch=compute_50,code=sm_50
    -gencode=arch=compute_52,code=sm_52
    -gencode=arch=compute_60,code=sm_60
    -gencode=arch=compute_61,code=sm_61
    -gencode=arch=compute_70,code=sm_70
    -gencode=arch=compute_75,code=sm_75
    -gencode=arch=compute_75,code=compute_75
    -ftz=true -lineinfo
    -std=c++11)

set(CUDA_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/lib)
if (EXISTS ${CUDA_TOOLKIT_ROOT_DIR}/samples/common/inc)
  list(APPEND CUDA_INCLUDES ${CUDA_TOOLKIT_ROOT_DIR}/samples/common/inc)
endif()

cuda_add_library(pot STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/double_env.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/gpupot.gpu.cu OPTIONS ${CUDA_OPTIONS})
target_include_directories(pot PUBLIC ${CUDA_INCLUDES})

add_library(gpunb_bonsai STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/gpunb.bonsai.cpp)
target_include_directories(gpunb_bonsai PUBLIC
  ${CUDA_INCLUDES}
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/Bonsai/include)

cuda_add_library(gpunb_velocity
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/gpunb.velocity.cu OPTIONS ${CUDA_OPTIONS})
target_include_directories(gpunb_velocity PUBLIC ${CUDA_INCLUDES})

add_library(irr_tree STATIC ${CMAKE_CURRENT_SOURCE_DIR}/irrlib/gpuirr.tree.avx.cpp)
add_library(irr_avx STATIC ${CMAKE_CURRENT_SOURCE_DIR}/irrlib/gpuirr.avx.cpp)

set(SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/double_env.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/stopwatch.c
  ${CMAKE_CURRENT_SOURCE_DIR}/energy2.f
  ${CMAKE_CURRENT_SOURCE_DIR}/phicor.f
  ${CMAKE_CURRENT_SOURCE_DIR}/checkl2.f
  ${CMAKE_CURRENT_SOURCE_DIR}/swap.f
  ${CMAKE_CURRENT_SOURCE_DIR}/scale.f
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/wtime.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/cxvpred.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fpcorr2.f
  ${CMAKE_CURRENT_SOURCE_DIR}/repair.f
  ${CMAKE_CURRENT_SOURCE_DIR}/sweep2.f
  ${CMAKE_CURRENT_SOURCE_DIR}/ksres3.f
  ${CMAKE_CURRENT_SOURCE_DIR}/check3.f
  ${CMAKE_CURRENT_SOURCE_DIR}/bhplot.f
  ${CMAKE_CURRENT_SOURCE_DIR}/kspinit.f
  ${CMAKE_CURRENT_SOURCE_DIR}/kspreg.f
  )

SET(KS_SOURCES
  nbody6.f;assess.f;bindat.f;binev.f;binout.f;binpop.f;
  block.f;bodies.f;brake.f;brake2.f;brake3.f;bsetid.f;chaos0.f;
  chaos.f;check.f;checkl.f;chrect.f;clint.f;cloud.f;cloud0.f;
  cmbody.f;cmcorr.f;cmfreg.f;coal.f;comenv.f;
  corerd.f;cputim.f;data.f;decide.f;define.f;deform.f;degen.f;
  delay.f;dgcore.f;dtchck.f;eccmod.f;ecirc.f;edot.f;efac2.f;
  efac3.f;expel.f;energy.f;escape.f;events.f;evolve.f;expand.f;
  fbulge.f;fclose.f;fcloud.f;fcorr.f;fdisk.f;fhalo.f;ficorr.f;findj.f;
  findm.f;flyby.f;fnuc.f;forceir13.f;fpcorr.f;fpoly1.f;fpoly2.f;
  freeze.f;gcinit.f;gcint.f;giant.f;giant3.f;gntage.f;grrad.f;
  hcorr.f;hiarch.f;hicirc.f;hidat.f;higrow.f;himax.f;himax2.f;
  himod.f;hipop.f;hirect.f;histab.f;hivel.f;hmdot.f;hmdot2.f;
  hotsys.f;hrdiag.f;hrplot.f;hut.f;hut2.f;iblock.f;imf.f;imfbd.f;
  imf2.f;impact.f;induce.f;inext.f;input.f;insert.f;instar.f;
  intgrt.f;intide.f;jacobi.f;kepler.f;kick.f;kick2.f;kozai.f;ksapo.f;
  kscorr.f;ksin2.f;ksmod.f;ksperi.f;
  ksrect.f;ksres.f;ksres2.f;
  kstide.f;lagr.f;lagr2.f;levels.f;magbrk.f;matrix.f;merge.f;
  merge2.f;mix.f;mloss.f;mlwind.f;modify.f;mrenv.f;mtrace.f;
  nbhist.f;nbint.f;nbpot.f;nbrem.f;nbrest.f;nbsort.f;nbtide.f;
  newtev.f;nstab.f;ntint.f;offset.f;orbit.f;peri.f;permit.f;
  pfac.f;plpot1.f;poti.f;premsf.f;proto_star.f;qtides.f;ran2.f;reflct.f;
  regint.f;remove.f;rename.f;reset.f;reset2.f;rkint.f;rl.f;roche.f;
  rpmax.f;rpmax2.f;rpmin.f;scale.f;search.f;setup.f;setup2.f;shrink.f;
  sort1.f;spiral.f;stability.f;star.f;start.f;stepk.f;steps.f;stumpf.f;
  subint.f;subsys.f;swap.f;sweep.f;synch.f;tail0.f;tcirc.f;tides.f;tides2.f;
  tides3.f;touch.f;trdot.f;trdot2.f;trflow.f;tstab.f;tstep.f;
  units.f;unpert.f;update.f;verify.f;xtrnl0.f;xtrnld.f;xtrnlf.f;
  xtrnlp.f;xtrnlt.f;xtrnlv.f;zare.f;zcnsts.f;zero.f;zfuncs.f;
  triple.f;derqp3.f;difsy3.f;erel3.f;extend.f;qpmod3.f;stabl3.f;
  stablz.f;start3.f;tperi.f;trans3.f;
  quad.f;derqp4.f;difsy4.f;endreg.f;erel4.f;ichain.f;newreg.f;newsys.f;
  qpmod4.f;rchain.f;rsort.f;stabl4.f;start4.f;status.f;trans4.f;
  brake4.f;ksrot.f;expand2.f;deform2.f;circ.f;spinup.f;qstab.f;
  newton.f;pnpert2.f
  )

if (USE_KSBLOCK)
  SET(KS_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../Block)
  LIST(APPEND KS_SOURCES ksintp.f;kspert2.f;kstrans2.f;stepi.f)
else()
  SET(KS_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode)
  LIST(APPEND KS_SOURCES kstran.f)
endif()

foreach(KS_SOURCE ${KS_SOURCES})
  find_file(SOURCE_FILE ${KS_SOURCE}
    PATHS ${CMAKE_CURRENT_SOURCE_DIR} ${KS_PATH}
    NO_DEFAULT_PATH)
  IF (NOT SOURCE_FILE)
    MESSAGE(FATAL_ERROR "${KS_SOURCE} not found")
  ENDIF()
  LIST(APPEND SOURCES ${SOURCE_FILE})
  unset(SOURCE_FILE CACHE)
endforeach()

SET(CHAIN_SOURCES
  chain.f;chstab.f;const.f;cstab2.f;cstab3.f;cstab4.f;
  erel.f;hpsort.f;inclin.f; ksphys.f;physks.f;r2sort.f;
  absorb.f;chdata.f;chfind.f;chinit.f;chmod.f;
  chpert.f;chpot.f;expel2.f;ghost.f;chterm.f;chterm2.f;
  reduce.f;reinit.f;renew.f;setsys.f;tchain.f;xtpert.f)

if (USE_ARCHAIN)
  SET(CHAIN_PATHS
    ${CMAKE_CURRENT_SOURCE_DIR}/../ARchain
    ${CMAKE_CURRENT_SOURCE_DIR}/../ARint)


  LIST(APPEND CHAIN_SOURCES
    arc.f;bhstab.f;bhstat.f;decorr.f;
    emax1.f;ffdot.f;grbin.f;infall.f;inject.f;pnpert.f;xtpert2.f)
else()
  SET(CHAIN_PATHS
    ${CMAKE_CURRENT_SOURCE_DIR}/../Chain
    ${CMAKE_CURRENT_SOURCE_DIR}/../Nchain)

  LIST(APPEND CHAIN_SOURCES
    cfuncs.f;chaos2.f;cstab5.f;derqp.f;difsy1.f;giant2.f;invert.f;
    qforce.f;qpmod.f;recoil.f;redraw.f;select.f;slow.f;stablc.f;
    swcond.f;switch.f;transk.f;transq.f;transx.f;vector.f;xtf.f;
    xtrnlu.f;ycopy.f;ysave.f)
endif()

foreach(CH_SOURCE ${CHAIN_SOURCES})
  find_file(SOURCE_FILE ${CH_SOURCE}
    PATHS ${CMAKE_CURRENT_SOURCE_DIR} ${CHAIN_PATHS}
    NO_DEFAULT_PATH)
  IF (NOT SOURCE_FILE)
    MESSAGE(FATAL_ERROR "${CH_SOURCE} not found")
  ENDIF()
  LIST(APPEND SOURCES "${SOURCE_FILE}")
  unset(SOURCE_FILE CACHE)
endforeach()

add_library(nbody STATIC ${SOURCES})

set(TREE_CHAIN_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/chfirr.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/kcpert.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/fchain.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/chlist.tree.f

#  ${CMAKE_CURRENT_SOURCE_DIR}/../ARint/chfirr.f
#  ${CMAKE_CURRENT_SOURCE_DIR}/../ARint/kcpert.f
#  ${CMAKE_CURRENT_SOURCE_DIR}/../ARint/fchain.f

  ${CMAKE_CURRENT_SOURCE_DIR}/xcpred.tree.f
  )
if (USE_ARCHAIN)
  LIST(APPEND TREE_CHAIN_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/cpertj.tree.f
    ${CMAKE_CURRENT_SOURCE_DIR}/cpertx.tree.f
    ${CMAKE_CURRENT_SOURCE_DIR}/chterm2.tree.f

#    ${CMAKE_CURRENT_SOURCE_DIR}/../ARint/cpertj.f
#    ${CMAKE_CURRENT_SOURCE_DIR}/../ARint/cpertx.f
    )
endif()

SET(GPU_CH_FILES
  chfirr.f;kcpert.f;fchain.f;chlist.f;xcpred.f;cpertj.f;cpertx.f;
  )
SET(GPU_CH_SOURCES)
foreach(CH_SOURCE ${GPU_CH_FILES})
  find_file(SOURCE_FILE ${CH_SOURCE}
    PATHS ${CMAKE_CURRENT_SOURCE_DIR} ${CHAIN_PATHS}
    NO_DEFAULT_PATH)
  IF (NOT SOURCE_FILE)
    MESSAGE(FATAL_ERROR "${CH_SOURCE} not found")
  ENDIF()
  LIST(APPEND GPU_CH_SOURCES "${SOURCE_FILE}")
  unset(SOURCE_FILE CACHE)
endforeach()

add_subdirectory(lib/Bonsai)
add_executable(nbody7.tree
  ${CMAKE_CURRENT_SOURCE_DIR}/start.f
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/knn.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/core.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/mydump.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/fpoly0.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/intgrt.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/adjust.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/steps.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/jpred.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/nbintp.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/nbint.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/tstep.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/search.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/output.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/rkick.f
  ${CMAKE_CURRENT_SOURCE_DIR}/gpuint.f
  ${CMAKE_CURRENT_SOURCE_DIR}/gpuint2.f
  ${CMAKE_CURRENT_SOURCE_DIR}/ksinit.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/ksreg.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/ksterm.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/kslist.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/cutoff.f
  ${CMAKE_CURRENT_SOURCE_DIR}/cutoff0.f
  ${CMAKE_CURRENT_SOURCE_DIR}/nblist.tree.f

  ${TREE_CHAIN_SOURCES}
  #${GPU_CH_SOURCES}
  ${CMAKE_CURRENT_SOURCE_DIR}/xcpred.tree.f

  ${CMAKE_CURRENT_SOURCE_DIR}/tpert.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/ksint.tree.f

  ${CMAKE_CURRENT_SOURCE_DIR}/jpred2.tree.f
  #${CMAKE_CURRENT_SOURCE_DIR}/jpred2.f

  #${CMAKE_CURRENT_SOURCE_DIR}/fpert.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/fpert.f

  #${CMAKE_CURRENT_SOURCE_DIR}/kspoly.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/kspoly.f

  ${CMAKE_CURRENT_SOURCE_DIR}/xvpred.tree.f
  #${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/xvpred.f

  ${CMAKE_CURRENT_SOURCE_DIR}/resolv.tree.f
  #${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/resolv.f

  #${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/kspred.f
  ${CMAKE_CURRENT_SOURCE_DIR}/kspred.tree.f

  ${CMAKE_CURRENT_SOURCE_DIR}/cmfirr.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/cmfirr2.tree.f
  #${CMAKE_CURRENT_SOURCE_DIR}/cmfirr.f
  #${CMAKE_CURRENT_SOURCE_DIR}/cmfirr2.f
  ${CMAKE_CURRENT_SOURCE_DIR}/kspert.tree.f
  #${CMAKE_CURRENT_SOURCE_DIR}/kspert.f
  #${CMAKE_CURRENT_SOURCE_DIR}/lib/cnbint.tree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/cnbint.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/auswert.f
  ${CMAKE_CURRENT_SOURCE_DIR}/auswert_iso.f
  ${CMAKE_CURRENT_SOURCE_DIR}/sort2.f
  ${CMAKE_CURRENT_SOURCE_DIR}/binlist.f
  ${CMAKE_CURRENT_SOURCE_DIR}/dcen.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/mdot.f
  )
add_executable(nbody7.leapfrog
  ${CMAKE_CURRENT_SOURCE_DIR}/fpoly0.leapfrog.f
  ${CMAKE_CURRENT_SOURCE_DIR}/start.leapfrog.f
  ${CMAKE_CURRENT_SOURCE_DIR}/intgrt.leapfrog.f
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/knn.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/core.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/mydump.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/adjust.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/steps.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/jpred.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/tstep.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/search.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/output.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/ksinit.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/ksreg.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/ksterm.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/kslist.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/cutoff.f
  ${CMAKE_CURRENT_SOURCE_DIR}/cutoff0.f
  ${CMAKE_CURRENT_SOURCE_DIR}/nblist.tree.f

  ${TREE_CHAIN_SOURCES}
  ${CMAKE_CURRENT_SOURCE_DIR}/xcpred.tree.f

  ${CMAKE_CURRENT_SOURCE_DIR}/tpert.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/ksint.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/jpred2.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/fpert.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/kspoly.f
  ${CMAKE_CURRENT_SOURCE_DIR}/xvpred.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/resolv.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/kspred.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/cmfirr.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/cmfirr2.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/kspert.tree.f
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/cnbint.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/auswert.f
  ${CMAKE_CURRENT_SOURCE_DIR}/auswert_iso.f
  ${CMAKE_CURRENT_SOURCE_DIR}/sort2.f
  ${CMAKE_CURRENT_SOURCE_DIR}/binlist.f
  ${CMAKE_CURRENT_SOURCE_DIR}/dcen.f
  ${CMAKE_CURRENT_SOURCE_DIR}/mdot.leapfrog.f
  )

add_executable(nbody7.gpu
  ${CMAKE_CURRENT_SOURCE_DIR}/start.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/core.f
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/cnbint.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/mydump.f
  ${CMAKE_CURRENT_SOURCE_DIR}/fpoly0.f
  ${CMAKE_CURRENT_SOURCE_DIR}/intgrt.omp.f
  ${CMAKE_CURRENT_SOURCE_DIR}/adjust.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/steps.f
  ${CMAKE_CURRENT_SOURCE_DIR}/jpred.f
  ${CMAKE_CURRENT_SOURCE_DIR}/nbintp.f
  ${CMAKE_CURRENT_SOURCE_DIR}/nbint.f
  ${CMAKE_CURRENT_SOURCE_DIR}/gpucor.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/tpert.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/tstep.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/search.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/output.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/ksterm.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/ksinit.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/ksreg.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/kslist.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/kspoly.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/ksint.f

  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/xvpred.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/kspred.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/fpert.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/resolv.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/nblist.f
  ${CMAKE_CURRENT_SOURCE_DIR}/cmfirr.f
  ${CMAKE_CURRENT_SOURCE_DIR}/cmfirr2.f
  ${CMAKE_CURRENT_SOURCE_DIR}/jpred2.f
  ${CMAKE_CURRENT_SOURCE_DIR}/kspert.f
  ${CMAKE_CURRENT_SOURCE_DIR}/../Ncode/mdot.f
  ${GPU_CH_SOURCES}
  )


target_link_libraries(nbody7.tree nbody gpunb_bonsai pot irr_tree bonsai)
target_link_libraries(nbody7.leapfrog nbody gpunb_bonsai pot irr_tree bonsai)
target_link_libraries(nbody7.gpu nbody gpunb_velocity pot irr_avx)

add_executable(getproj ${CMAKE_CURRENT_SOURCE_DIR}/lib/getproj.c)
target_link_libraries(getproj m)
