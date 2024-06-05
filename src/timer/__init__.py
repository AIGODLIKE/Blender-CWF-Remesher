import bpy
import traceback
from queue import Queue
from typing import Any


class Timer:
    TimerQueue = Queue()

    @staticmethod
    def put(delegate: Any):
        Timer.TimerQueue.put(delegate)

    @staticmethod
    def executor(t):
        if type(t) in {list, tuple}:
            t[0](*t[1:])
        else:
            t()

    @staticmethod
    def run1():
        return Timer.run_ex(Timer.TimerQueue)

    @staticmethod
    def run_ex(queue: Queue):
        while not queue.empty():
            t = queue.get()
            # Timer.executor(t)
            try:
                Timer.executor(t)
            except Exception:
                traceback.print_exc()
            except KeyboardInterrupt:
                ...
        return 0.016666666666666666

    @staticmethod
    def clear():
        while not Timer.TimerQueue.empty():
            Timer.TimerQueue.get()

    @staticmethod
    def wait_run(func):
        def wrap(*args, **kwargs):
            q = Queue()

            def wrap_job(q):
                try:
                    res = func(*args, **kwargs)
                    q.put(res)
                except Exception as e:
                    q.put(e)

            Timer.put((wrap_job, q))
            res = q.get()
            if isinstance(res, Exception):
                raise res
            return res

        return wrap

    @staticmethod
    def reg():
        bpy.app.timers.register(Timer.run1, persistent=True)

    @staticmethod
    def unreg():
        Timer.clear()
        try:
            bpy.app.timers.unregister(Timer.run1)
        except Exception:
            ...


def register():
    Timer.reg()


def unregister():
    Timer.unreg()
