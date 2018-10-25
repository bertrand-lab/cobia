from setuptools import setup
import cobia

setup(name = 'cobia',
      version = cobia.__version__,
      description = 'cobia: predicting cofragmentation induced bias',
      author = 'Scott McCain',
      url = 'https://github.com/bertrand-lab/cobia',
      author_email = "j.scott.mccain@dal.ca",
      license = 'GPL3',
      packages = ['cobia'],
      entry_points = {'console_scripts': ['cobia=cobia.cobia:main']},
      zip_safe = False)
