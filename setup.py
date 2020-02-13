from setuptools import setup, Extension
import sys


class _deferred_pybind11_include(object):
    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


include_dirs = ['./', _deferred_pybind11_include(True),
                _deferred_pybind11_include()]

extra_compile_args = ['--std=c++11', '-march=native', '-O3', '-fvisibility=hidden', '-fopenmp']
python_module_link_args = extra_compile_args
define_macros = []


def get_extension_modules():
    return [Extension('sph_tools',
                      language='c++',
                      sources=['sph_tools/module.cpp', 'sph_tools/h_proj_3d_core.c', 'sph_tools/kernel.c'],
                      depends=['sph_tools/h_proj.h'],
                      include_dirs=include_dirs,
                      define_macros=define_macros,
                      extra_compile_args=extra_compile_args,
                      extra_link_args=python_module_link_args)]


setup(name='sph_tools',
      version='0.0.1',
      description='Basic tools for plotting and analysing SPH data',
      include_package_data=True,
      author='Peter Bell',
      author_email='peterbelll10@live.co.uk',
      packages=[],
      setup_requires=['numpy>=1.15.0', 'pybind11>=2.2.4'],
      ext_modules=get_extension_modules(),
      install_requires=['numpy>=1.15.0', 'pybind11>=2.2.4']
      )
