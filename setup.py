from setuptools import setup, Extension
from Cython.Build import cythonize

# Define the Cython extensions
extensions = [
    Extension(
        "sequence_clustering.dsu",
        ["src/sequence_clustering/dsu.pyx"],
        include_dirs=[],
        language="c",
    ),
    Extension(
        "sequence_clustering.utils",
        ["src/sequence_clustering/utils.pyx"],
        include_dirs=[],
        language="c",
    )
]

setup(
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            'language_level': 3,
            'boundscheck': True,  # Enable bounds checking for safety
            'wraparound': False,
            'cdivision': True,
            'nonecheck': False,
        }
    ),
    zip_safe=False,
)
