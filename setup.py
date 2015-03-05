try:
    from setuptools import setup, Command
except ImportError:
    from distutils.core import setup, Command

import pkg_resources

##
# Simple command for running tests.
#
class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import sys,subprocess
        errno = subprocess.call([sys.executable, 'runtests.py'])
        raise SystemExit(errno)

setup(
    name="svsim",
    version="0.2dev",
    description = "Tool for simelating structural variations",
    long_description = open( "README.md", "r" ).read( ),
    install_requires = [
        'click', 
        'pyfasta',
        'numpy',
        'pysam'
    ],
    package_dir = { 
        "svsim" : "svsim" 
    },
    packages = [
        'svsim', 
        'svsim.reads',
        'svsim.commands'
    ],
    package_data= { 
        'svsim': [ 'data/*.mconf' ] 
    },
    scripts = [
        'scripts/svsim'
    ],
    license="Modified BSD",
    url="https://www.github.com/fadern/svsim",
    maintainer="The svest team.",
    maintainer_email="mattias.franberg@scilifelab.se",
    cmdclass = { 
        "test" : PyTest 
    }
)
