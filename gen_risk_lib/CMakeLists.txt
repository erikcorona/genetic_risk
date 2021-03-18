project(gen_risk_lib)

set(HEADER_FILES GWAS.hxx FlatFile.hxx)

add_library(gen_risk_lib STATIC ${HEADER_FILES})