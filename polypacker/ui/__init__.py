from .headless import CLI
from .graphical import GUI


def uipicker(uitype):
    if uitype.lower() in {'graphical', 'gui'}:
        return GUI
    elif uitype.lower() in {'headless', 'cli'}:
        return CLI
    else:
        return None
