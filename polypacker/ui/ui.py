from abc import ABC, abstractmethod


class UI(ABC):

    # ============================== #
    #              INIT              #
    # ============================== #

    def __init__(self, ITER_max, ITER_upd=1):
        self.ITER_max = ITER_max  # positive integer
        self.ITER_upd = ITER_upd  # positive integer, cannot be zero (set to inf instead!)

    @abstractmethod
    def init(self):
        pass

    # ============================== #
    #              STOP              #
    # ============================== #

    @property
    @abstractmethod
    def STOP(self):
        """UI_Action to stop."""
        return None  # default None = always False

    @property
    def was_stopped(self):
        if self.STOP is None:
            return False
        else:
            return self.STOP.is_active

    # ============================== #
    #             UPDATE             #
    # ============================== #

    @property
    @abstractmethod
    def UPDATE(self):
        """UI_Action to update."""
        return None  # default None = always False

    def needs_update(self, ITER, **kwargs):
        if self.UPDATE is None or self.was_stopped:
            UPDATE = False
        else:
            UPDATE = self.UPDATE.is_active
        return UPDATE or update_needed(ITER, self.ITER_max, every_n=self.ITER_upd, **kwargs)

    def do_update(self, ITER, polygons=None, density=None, collisions=None, **kwargs):
        if not self.needs_update(ITER):  # double-check
            return  # terminate early
        self._do_update(ITER, polygons, density, collisions, **kwargs)  # pass through

    @abstractmethod
    def _do_update(self, ITER, polygons, density, collisions, **kwargs):
        pass


class UI_Action(ABC):

    @property
    @abstractmethod
    def is_active(self):
        """Boolean indicating if the Action is active."""
        return False

    @abstractmethod
    def reset(self):
        """Reset the Action."""
        pass


def update_needed(n, n_max, every_n=10, at_start=True, at_end=True):
    """Check if UI needs updating based on iteration number."""
    is_start = n==0 if at_start else False  # always at the beginning
    is_end = n==n_max-1 if at_end else False  # always at the end
    is_n = (n%every_n)==0  # every n
    return is_start or is_n or is_end
