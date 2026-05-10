from setuptools import setup, Extension, find_packages
import os
import glob

sources = []
sources += [f for f in glob.glob("src/*.cpp") if os.path.basename(f) != "factorizer.cpp"]
sources += glob.glob("src/*.pyx")

root_dir = os.path.abspath(os.path.dirname(__file__))
    
ext = Extension("factorizer", 
    sources = sources,
    language = "c++",
    extra_compile_args = ["-v", "-std=c++14", "-Wall", "-O3", "-lboost_system"],
    extra_link_args = ["-std=c++14"]
)

with open(os.path.join(root_dir, 'README.md'), "r") as fp:
    long_description = fp.read()
with open(os.path.join(root_dir, 'requirements.txt'), "r") as fp:
    install_requires = fp.read().splitlines()

setup(
    name = "factorizer",
    version = "0.9.6",
    author = "Fulltea",
    author_email = "rikuta@furutan.com",
    long_description = long_description,
    long_description_content_type="text/markdown",
    url = "https://github.com/FullteaOfEEIC/factorizer",
    packages = find_packages(where="src"),
    package_dir = {
        "factorizer": "src"
    },
    install_requires = install_requires,
    ext_modules = [ext]
)

if os.path.exists(os.path.join("src", "factorizer.cpp")):
    os.remove(os.path.join("src", "factorizer.cpp"))
