from conan import ConanFile
from conan.tools.layout import basic_layout
from conan.tools.files import copy

class EZArpackConan(ConanFile):
    name = "ezarpack"
    version = "1.1"
    package_type = "header-library"

    # Metadata
    license = "MPL-2.0"
    author = "Igor Krivenko (iskrivenko@proton.me)"
    url = "https://github.com/krivenko/ezARPACK.git"
    homepage = "https://krivenko.github.io/ezARPACK/"
    description = "A C++11 wrapper around ARPACK-NG compatible with Eigen3, Armadillo, Blaze, xtensor and other matrix algebra libraries"
    topics = ("arpack", "eigen", "boost", "armadillo", "blaze", "xtensor")

    exports_sources = "include/*"

    implements = ["auto_header_only"]

    def layout(self):
        basic_layout(self)

    def package(self):
        copy(self, "include/*", self.source_folder, self.package_folder)

    def package_info(self):
        self.cpp_info.system_libs = ["arpack", "lapack", "blas"]
        self.cpp_info.bindirs = []
        self.cpp_info.libdirs = []


