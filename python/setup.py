from distutils.core import setup
from Cython.Distutils import build_ext
from distutils.extension import Extension

sources_list = [
    # "ising_impl.pxd",
    "ising_impl.pyx",
    "../source/ising_impl.cpp",
]

setup(
    ext_modules=[
        Extension(
            "ising_sim",
            sources=sources_list,
            language="c++",
            extra_compile_args=["-std=c++11"],
            cython_directives={"language_level": "3"},
            include_dirs=["../include", "../source"],
        )
    ],
    cmdclass={"build_ext": build_ext},
)
