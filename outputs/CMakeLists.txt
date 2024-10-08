CMAKE_MINIMUM_REQUIRED(VERSION 3.10)
project(outputs)

enable_language(C)
add_compile_options(-std=c11 -pedantic -Wall -Wextra -Werror -pedantic-errors)
find_package(OpenMP)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(outputs
  ${CMAKE_CURRENT_SOURCE_DIR}/src/bound.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/clmass.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/data.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/discard.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dot.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/forceir13.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/fpstart.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/freefrm.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/frmhdr.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/frmlst.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/frmmeta.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/frmsz.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/hdr.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mag2.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/meta.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/numfrm.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/oerrno.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/oerrstr.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/oerror.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ofail.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/opnout.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/rdfrm.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/reclen.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/readrec.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/rtide.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/skpfrm.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/vsub.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/vmul.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/cross.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/galforce.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/forcebv15.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/center.c
  )
target_link_libraries(outputs m)

add_executable(count_frames
  ${CMAKE_CURRENT_SOURCE_DIR}/count_frames.c)
target_link_libraries(count_frames outputs)


add_executable(list_frames
  ${CMAKE_CURRENT_SOURCE_DIR}/list_frames.c)
target_link_libraries(list_frames outputs)

add_executable(rtide
  ${CMAKE_CURRENT_SOURCE_DIR}/rtide.c)
target_link_libraries(rtide outputs)

add_executable(boundm
  ${CMAKE_CURRENT_SOURCE_DIR}/boundm.c)
target_link_libraries(boundm outputs)

add_executable(info
  ${CMAKE_CURRENT_SOURCE_DIR}/info.c)
target_link_libraries(info outputs)

add_executable(orbit
  ${CMAKE_CURRENT_SOURCE_DIR}/orbit.c)
target_link_libraries(orbit outputs)

add_executable(run_info
  ${CMAKE_CURRENT_SOURCE_DIR}/run_info.c)
target_link_libraries(run_info outputs)

add_executable(snapshot
  ${CMAKE_CURRENT_SOURCE_DIR}/snapshot.c)
target_link_libraries(snapshot outputs)

add_executable(stream
  ${CMAKE_CURRENT_SOURCE_DIR}/stream.c)
target_link_libraries(stream outputs)

find_package(OpenMP)
add_executable(vel
  ${CMAKE_CURRENT_SOURCE_DIR}/vel.c)
target_link_libraries(vel outputs OpenMP::OpenMP_C)
