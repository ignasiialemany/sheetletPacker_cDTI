import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

from .ui import UI, UI_Action, update_needed
from ..visualize import plot_polygon, update_patch

class GUI(UI):

    # ============================== #
    #              INIT              #
    # ============================== #

    def init(self, polygons, region=None, *, polygonscolor='tab:blue', regioncolor='tab:red', **kwargs):

        # create the figure and subplots
        fig, axs = plt.subplots(ncols=3, nrows=2)
        fig.canvas.set_window_title(kwargs.pop('window_title', ''))

        # create one large subplot from smaller ones
        grid = axs[0, 0].get_gridspec()  # get grid from top left reference
        [axs[row, col].remove() for row in range(2) for col in range(2)]
        axes = dict(poly=fig.add_subplot(grid[:2, :2]),
                    dens=axs[0, -1],
                    coll=axs[1, -1])

        # show all polygons
        axes['poly'].set_title('polygons')
        unit = kwargs.pop('distance_unit', '-')
        axes['poly'].set_xlabel('x ({})'.format(unit))
        axes['poly'].set_ylabel('y ({})'.format(unit))
        axes['poly'].set_aspect('equal', adjustable='datalim')

        # density for convergence
        axes['dens'].set_title('density')
        axes['dens'].set_xlabel('iteration #')
        axes['dens'].set_ylabel('density (%)')
        axes['dens'].set_xlim([0, max(self.ITER_max, 2)-1])
        axes['dens'].set_ylim([0, 1])

        # number of collisions
        axes['coll'].set_title('collisions')
        axes['coll'].set_xlabel('iteration #')
        axes['coll'].set_ylabel('# of collisions')
        axes['coll'].set_xlim([0, max(self.ITER_max, 2)-1])
        axes['coll'].autoscale(axis='y')

        # plot region and initial polygon data
        def plot_poly(poly, color, axis):
            face = poly.color if hasattr(poly, 'color') else color
            return plot_polygon(poly, axis=axis, facecolor=face)
        ax = axes['poly']
        self._region = plot_poly(region, regioncolor, ax) if region is not None else None
        self._patches = [plot_poly(p, polygonscolor, ax) for p in polygons]

        # plot line data (hidden due to nan value)
        empty_nan = np.full(self.ITER_max, np.nan)
        # need to .copy() the data, otherwise plot data is linked
        self._density = axes['dens'].plot(empty_nan.copy(), 'r.-')[0]
        self._collisions = axes['coll'].plot(empty_nan.copy(), 'b.-')[0]

        # create axes for buttons at [left, bottom, width, height]
        h, w = 0.05, 0.09
        axes['stop'] = plt.axes([0.01, 0.01, w, h])
        axes['update'] = plt.axes([0.11, 0.01, w, h])
        self._STOP = GUI_ToggleButton(axes['stop'], 'Stop')
        self._UPDATE = GUI_ToggleButton(axes['update'], 'Update')

        # show
        fig.tight_layout()
        plt.ion()  # turn on interactive mode (shows)
        plt.show()  # does not lock thanks to ion()
        self._fig = fig
        self._axes = axes

        return axes  # return the axes dict in case user wants to modify

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
        return self._UPDATE

    def _do_update(self, ITER, polygons=None, density=None, collisions=None, *, wait_time=0.001):

        # polygons
        if polygons is not None:
            for patch, poly in zip(self._patches, polygons):
                update_patch(patch, poly)

        # density & collisions (latter needs to update the limits dynamically)
        def update_line(i, line, new, relim_ax=None):
            if new is not None:
                ydat = line.get_ydata()
                ydat[i] = new
                line.set_ydata(ydat)
                if relim_ax is not None:
                    relim_ax.relim()
                    relim_ax.autoscale_view()
        update_line(ITER, self._density, density)
        update_line(ITER, self._collisions, collisions, relim_ax=self._axes['coll'])

        # process
        self._fig.canvas.draw()
        self._fig.canvas.flush_events()
        plt.draw()
        plt.pause(wait_time)

        # done
        self.UPDATE.reset()


class GUI_ToggleButton(UI_Action):
    """ToggleButton action class for GUI.
    See https://stackoverflow.com/a/28760906
    """

    def __init__(self, axes, label):
        self._toggled = False  # default value
        self._button = Button(axes, label)
        self._button.on_clicked(self._pressed)

    def _pressed(self, event):
        self._toggled = True

    @property
    def is_active(self):
        return self._toggled

    def reset(self):
        self._toggled = False
