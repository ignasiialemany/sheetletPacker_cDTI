import os.path
import math

from .ui import UI, UI_Action, update_needed


class CLI(UI):

    def init(self, *args, **kwargs):
        self._STOP = CLI_StopAction()

    # ============================== #
    #              STOP              #
    # ============================== #

    @property
    def STOP(self):
        return self._STOP

    # ============================== #
    #             UPDATE             #
    # ============================== #

    @property
    def UPDATE(self):
        return None

    def _do_update(self, ITER, polygons=None, density=None, collisions=None, *, header_every=20, fillchar='-'):
        # progress
        L = math.floor(math.log10(self.ITER_max))+1  # number of digits for ITER
        fmt = 'ITER: %'+str(L)+'d of %'+str(L)+'d = %5.1f%%'
        dat = (ITER+1, self.ITER_max, (ITER+1)/self.ITER_max*100)
        nchar = 6+L+4+L+3+5+1  # number of characters in fmt. Do NOT use len(fmt)
        hdr = '#' + 'Progress'.center(nchar-1, fillchar)
        # extra
        sep = ' | '  # separate from previous text
        w = 20
        # density
        hdr += sep.replace(' ', fillchar) + 'density'.center(w, fillchar)
        fmt += sep + ('%'+str(w)+'g' if density else '%s')
        dat += (density if density else ' '*w,)
        # collisions
        hdr += sep.replace(' ', fillchar) + 'num_hits'.center(w, fillchar)
        fmt += sep + ('%'+str(w)+'d' if collisions else '%s')
        dat += (collisions if collisions else ' '*w,)
        # print
        fmt += ' |'
        hdr += fillchar + '|'
        if update_needed(ITER//self.ITER_upd, self.ITER_max//self.ITER_upd, every_n=header_every, at_end=False):
            print(hdr)
        print(fmt % dat)


class CLI_StopAction(UI_Action):

    @property
    def is_active(self):
        return os.path.isfile('STOP') or os.path.isfile('ABORT')

    def reset(self):
        if os.path.isfile('STOP'):
            os.remove('STOP')
        if os.path.isfile('ABORT'):
            os.remove('ABORT')
