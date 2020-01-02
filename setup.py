from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='AREAMpy',
      version='0.1',
      description='Amercian Railway Engineering and Maintenance-of-Way Association in Python',
      long_description=readme(),
      license='GPL 3.0',
      author='mlw',
      url='https://github.com/mwhit74/AREMApy',
      packages=[],
      install_requires=[])
