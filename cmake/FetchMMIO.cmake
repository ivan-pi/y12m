# cmake/FetchMMIO.cmake
#
# Downloads mmio.f from NIST into the build tree.
# Sets MMIO_F to the full path on success, or raises a FATAL_ERROR.
#
# Update MMIO_SHA256 if NIST ever updates the file:
#   curl -fsSL https://math.nist.gov/MatrixMarket/mmio/f/mmio.f | sha256sum

set(MMIO_URL    "https://math.nist.gov/MatrixMarket/mmio/f/mmio.f")
set(MMIO_SHA256 "f3146809daffd47587ecd0fa8af0187354a34bb41d3740cd71aab33e505d99a1")
set(MMIO_F      "${CMAKE_BINARY_DIR}/mmio/mmio.f")

file(DOWNLOAD "${MMIO_URL}" "${MMIO_F}"
    EXPECTED_HASH SHA256=${MMIO_SHA256}
    TLS_VERIFY ON
    STATUS _status
)
list(GET _status 0 _code)
if(NOT _code EQUAL 0)
    list(GET _status 1 _error)
    message(FATAL_ERROR "Failed to download mmio.f: ${_error}")
endif()