from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='rcssim',
    version='1.0.0',
    author='Tanner Dixon',
    author_email='tcdixon44@gmail.com',
    description='Functions for simulating the signal processing and logic operations of the Medtronic RC+S',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/Weill-Neurohub-OPTiMaL/rcs-simulation',
    project_urls = {
        "Bug Tracker": "https://github.com/Weill-Neurohub-OPTiMaL/rcs-simulation/issues"
    },
    packages=find_packages(),
    install_requires=['numpy', 'pandas', 'scipy', 'matplotlib'],
)