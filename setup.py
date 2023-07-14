from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension("orfclust.classes.treader", ["orfclust/classes/treader.pyx"]),
    Extension("orfclust.classes.transcript", ["orfclust/classes/transcript.pyx"]),
    Extension("orfclust.classes.txgroup", ["orfclust/classes/txgroup.pyx"]),
]

setup(
    name="ORFclust",
    version="0.0.1",
    author="Ales Varabyou, Zayn Zaidi",
    author_email="ales.varabyou@jhu.edu",
    description="Clustering transcriptome into translatome",
    url="https://github.com/alevar/orfclust",
    install_requires=["argparse","unittest"],
    python_requires='>=3.6',
    packages=['orfclust'],
    entry_points={'console_scripts': ['orfclust = orfclust.main:main'], },
    package_dir={"orfclust":"orfclust"},
    ext_modules=cythonize(extensions),
)
