from multiprocessing import Process, Queue
import os
import argparse
from vgg16_worker import Vgg16Worker


class Scheduler:
    def start(self, xfilelst):
        # put all of files into queue
        for xfile in xfilelst:
            self._queue.put(xfile)

        # add a None into queue to indicate the end of task
        self._queue.put(None)
        
#%%
class A(object):
    def __init__(self):
        self.x = 'Hello'

    def method_a(self, foo):
        print(self.x + ' ' + foo)
