import os


def neoart_sources(makefile):
    lines = open(makefile, 'r').readlines()

    files = []
    read = False
    for l in lines:
        if l.startswith('NEOART_SOURCES'): read = True
        if l.strip() == '': read = False
        if read: files.append(l.strip())

    files = [f.rstrip('\\') for f in files[1:]]

    path = os.path.dirname(makefile)
    files = [os.path.join(path, f) for f in files]
    return files


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('neoart', parent_package, top_path)

    neoart_files = neoart_sources('../src/Makefile-neoart.include')
    neoart_files += ['src/get_elem_config.f90', 'src/wrap_neoart.pyf']
    config.add_extension('_neoart', neoart_files)

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
