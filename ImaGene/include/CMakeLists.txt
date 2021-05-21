# Specifies how to install header files
file(GLOB_RECURSE HEADERS_H ImaGene *.h)
file(GLOB_RECURSE HEADERS_IH ImaGene *.ih)
list(LENGTH HEADERS_H NBH)
list(LENGTH HEADERS_IH NBIH)
message("Counted " ${NBH} " C++ headers .h")
message("Counted " ${NBIH} " C++ inline headers .ih")


