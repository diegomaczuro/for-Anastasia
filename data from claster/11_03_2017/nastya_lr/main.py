# -*- coding: utf_8 -*-
from postprocessing import *


def main():
    maps = Maps(FILE_MESH, 'out.vtu')
    maps._set_apd_treshold(-85.9654923732, 32.7753037234)
    maps.update()


if __name__ == '__main__':
    main()