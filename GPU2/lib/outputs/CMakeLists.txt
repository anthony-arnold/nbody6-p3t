CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

enable_language(C)

add_library(outputs
  ${CMAKE_CURRENT_SOURCE_DIR}/src/numfrm.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/oerrno.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/oerrstr.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/oerror.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ofail.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/opnout.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/reclen.c)


add_executable(count_frames
  ${CMAKE_CURRENT_SOURCE_DIR}/count_frames.c)
target_link_libraries(count_frames outputs)
