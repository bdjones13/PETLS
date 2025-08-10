from skbuild import setup

setup(
    name="petls",
    version="0.0.10",
    description="A Python Library to Compute Persistent Topological Laplacians",
    author="Benjamin Jones",
    license="",
    packages=["petls"],
    package_dir={"": "src"},
    cmake_install_dir="src/petls",
    python_requires=">=3.7",
)