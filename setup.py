from distutils.core import setup

setup(name='chesweet',
      version='0.0.1',
      description='Chemical shifts for glycans',
      author='Pablo Garay and Osvaldo Martin',
      author_email='garaypablo01@gmail.com, aloctavodia@gmail.com',
      url='https://github.com/BIOS-IMASL/chesweet',
      include_package_data = True,
      packages=['chesweet'],
      install_requires=['numpy', 'scipy'],
)
