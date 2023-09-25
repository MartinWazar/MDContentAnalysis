import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="MDContentAnalysis",
    version="1.0.4",
    author="Martin Wazar Eriksen",
    author_email="martin.wazar+GitHub@gmail.com",
    description="Package for analysing content of protein channels/tunnel in molecular dynamics simulations. Analyses trajectories. Built upon MDAnalysis, and uses Hole 2.0 or Mole 2.5",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MartinWazar/MDContentAnalysis",
    packages=setuptools.find_packages(),
    install_requires=['MDAnalysis'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)