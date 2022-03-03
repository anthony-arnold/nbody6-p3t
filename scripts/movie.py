#!/usr/bin/env python

import collections
import sys
import numpy
import warnings
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegFileWriter
from matplotlib.path import Path

ZMBAR = 0.626032
EX = Path([(-1, -1), (1, 1), (1, -1), (-1, 1)],
          [Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO])

EXTENTS = [-30, 30]
NORM = matplotlib.colors.Normalize(0, 10)

class Frame:
    def __init__(self, t):
        self.t = (t + ' ' * 11)[:11]
        self.dat = None

    def __repr__(self):
        l = ' empty'
        if self.dat is not None:
            l = ''
        return f'<Frame {self.t}{l}>'

def intensity(ax, frame, ext, bins, gal=True):
    # Get extents
    r = numpy.array([ext, ext])
    m = frame.dat[:,0]

    # Convert positions to kpc
    if not gal:
        x = frame.dat[:,4] / 1000
        y = frame.dat[:,5] / 1000

    else:
        x = frame.dat[:,1] / 1000
        y = frame.dat[:,2] / 1000


    h, _, _, i = ax.hist2d(x, y, bins=bins, weights=m, range=r,
                           cmap=matplotlib.cm.hot,
                           norm=NORM)

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', 'divide by zero encountered in log10')
        i.set_array(numpy.log10(h * 2500))
    return i

if __name__ == '__main__':
    px = 1/plt.rcParams['figure.dpi']

    plt.rcParams["figure.autolayout"] = True
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size'] = 18
    plt.style.use('dark_background')
    fig, ax = plt.subplots(ncols=2, figsize=(1280*px, 720*px))

    ax[0].set(xlim=EXTENTS, ylim=EXTENTS)
    ax[0].set_xlabel('x / kpc')
    ax[0].set_ylabel('y / kpc')

    cen = ax[0].scatter([0], [0], s=15, c='#ffffff', marker=EX)
    txt = ax[0].text(-29.9, 30.5, ' ' * 11)

    ax[1].set_xlabel('x / kpc')
    ax[1].set_ylabel('y / kpc')
    meshes = []
    cbs = []

    def get_frames():
        frames = collections.deque([
            Frame(sys.stdin.readline().rstrip())
        ])
        def snap():
            for line in sys.stdin:
                if line.startswith('T'):
                    frames.append(Frame(line.rstrip()))
                    return

                yield line

        while len(frames) > 0:
            frame = frames.popleft()
            frame.dat = numpy.loadtxt(snap())
            yield frame

    def render_frame(frame):
        print(frame)
        txt.set_text(f'{frame.t} Myrs')

        for mesh in meshes:
            mesh.remove()

        meshes[:] = [
            intensity(ax[0], frame, bins=3000, ext=EXTENTS),
            intensity(ax[1], frame, bins=100, ext=[-1, 1], gal=False)
        ]

        ax[0].set(xlim=EXTENTS, ylim=EXTENTS)
        ax[1].set(xlim=[-1,1], ylim=[-1,1])

        if len(cbs) == 0:
            cbs.append(
                fig.colorbar(ax=ax[1], mappable=None, cmap=matplotlib.cm.hot, norm=NORM,
                             label=r'$\log_{10}$ M$_{\odot}$ kpc$^{-2}$')
            )

    metadata = dict(title='NGC1904 Cartesian Projection')
    writer = FFMpegFileWriter(fps=24, metadata=metadata)
    with writer.saving(fig, 'movie.mp4', plt.rcParams['figure.dpi']):
        for frame in get_frames():
            render_frame(frame)
            writer.grab_frame()

        # Write the last frame out a few times.
        for _ in range(10):
            writer.grab_frame()
