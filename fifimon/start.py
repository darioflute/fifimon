#!/usr/bin/env python
from fifimon import mainwindow
import multiprocessing as mp

def main():
    mainwindow.main()

if __name__ == "__main__":
    mp.set_start_method('forkserver')
    main()    
