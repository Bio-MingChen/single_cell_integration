#!/usr/bin/env python
# -*- coding=utf-8 -*-

from memory_profiler import profile
from time import sleep

@profile
def test_profile():
    print('start')
    sleep(1)
    print('end')

@profile
def my_func():
    print('start')
    a = [1] * (10 ** 6)
    b = [2] * (2 * 10 ** 7)
    del b
    return a

if __name__ == "_main__":
    test_profile()
    my_func()