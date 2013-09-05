from distutils.core import setup, Command

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
    version="0.1dev",
    package_dir = { "svsim" : "svsim" },
    packages = ['svsim', 'svsim.reads'],
    package_data= { 'svsim': [ 'data/*.mconf' ] },
    license="Modified BSD",
    long_description = open( "README.md", "r" ).read( ),
    url="https://www.github.com/fadern/svsim",
    maintainer="The svest team.",
    maintainer_email="mattias.franberg@scilifelab.se",
    cmdclass = { "test" : PyTest }
)
