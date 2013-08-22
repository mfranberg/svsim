from distutils.core import setup

setup(
    name="svsim",
    version="0.1dev",
    package_dir = { "svsim" : "svsim" },
    packages = ['svsim', 'svsim.reads'],
    package_data= { 'svsim': [ 'data/*.mconf' ] },
    license="Modified BSD",
    long_description = open( "README.md", "r" ).read( ),
    url="http://www.example.com/",
    maintainer="The svest team.",
    maintainer_email="mattias.franberg@scilifelab.se"
)
