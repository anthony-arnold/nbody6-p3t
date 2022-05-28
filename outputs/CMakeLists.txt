CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(outputs)

enable_language(C)
add_compile_options(-std=c11 -pedantic -Wall -Wextra -Werror -pedantic-errors)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(outputs
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
  ${CMAKE_CURRENT_SOURCE_DIR}/src/vsub.c)
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
