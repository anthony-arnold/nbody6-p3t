#!/usr/bin/env python

"""
A script for turning NBODY6+P3T snapshots into a movie!

Expected input format is as follows:

    FILE: SNAPSHOT*
    SNAPSHOT: HEADER "\n" (PARTICLE "\n")*
    HEADER: "T = " FLOAT
    PARTICLE: FLOAT{7}


This format matches the streaming output of the `getproj` program:

    $ getproj <pos-file> stream | ./movie.py

The output of this script will be a 1280x720 MP4 @ 24fps.
"""

import collections
import sys
import numpy
import warnings
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegFileWriter
from matplotlib.path import Path

EX = Path([(-1, -1), (1, 1), (1, -1), (-1, 1)],
          [Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO])
NORM = matplotlib.colors.LogNorm(2, 9)

# TODO: Read these from input
ZMBAR = 0.626032
EXTENTS = [-30, 30]

class Frame:
    """
    A frame captures one snapshot from the input.
    """
    def __init__(self, t):
        """
        Initialise an empty frame with a specific timestamp.
        """
        self.t = (t + ' ' * 11)[:11]
        self.dat = None

    def __repr__(self):
        l = ' empty'
        if self.dat is not None:
            l = ''
        return f'<Frame {self.t}{l}>'

def intensity(ax, frame, ext, bins, gal=True):
    """
    Generate a 2D histogram for the given frame.
    The final weights will be log10(h) where h is the mass
    density per square kiloparsec in solar masses.

    :param ax: The Axes to draw on.
    :param frame: The Frame containing snapshot data.
    :param ext: The extents to sum, in the form [min, max].
    :param bins: The number of bins in each direction (see MPL hist2d).
    :param gal: If true, coordinates will be relative the the galaxy centre.
                otherwise will be relative to the cluster centre.
    :return: The Artist associated with the histogram.
    """
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
        i.set_array(numpy.log10(h * 2500)) # TODO: Parameterise this constant
    return i

if __name__ == '__main__':
    px = 1/plt.rcParams['figure.dpi']

    plt.rcParams["figure.autolayout"] = True
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size'] = 18
    plt.style.use('dark_background')

    fig, ax = plt.subplots(ncols=2, figsize=(1920*px, 1080*px))

    ax[0].set(xlim=EXTENTS, ylim=EXTENTS)
    ax[0].set_xlabel('x / kpc')
    ax[0].set_ylabel('y / kpc')

    # Draw an X at the galaxy centre.
    cen = ax[0].scatter([0], [0], s=15, c='#ffffff', marker=EX)

    # Write the time.
    txt = ax[0].text(-29.9, 30.5, ' ' * 11)

    ax[1].set_xlabel('x / kpc')
    ax[1].set_ylabel('y / kpc')
    meshes = []
    cbs = []

    def get_frames():
        """
        Stream the frames out one at a time.
        """
        frames = collections.deque([
            Frame(sys.stdin.readline().rstrip())
        ])
        def snap():
            """
            Read each line from a frame until the next frame header
            is encountered.
            """
            for line in sys.stdin:
                if line.startswith('T'):
                    frames.append(Frame(line.rstrip()))
                    return

                yield line

        # Read each frame from the input stream.
        while len(frames) > 0:
            frame = frames.popleft()
            frame.dat = numpy.loadtxt(snap())
            yield frame

    def render_frame(frame):
        print(frame)  # Report progress.
        txt.set_text(f'{frame.t} Myrs')

        # Remove old histograms from the plot.
        for mesh in meshes:
            mesh.remove()

        # Create new histograms.
        meshes[:] = [
            intensity(ax[0], frame, bins=3000, ext=EXTENTS),
            intensity(ax[1], frame, bins=100, ext=[-1, 1], gal=False)
        ]

        # Limits seem to need to be reset each time.
        ax[0].set(xlim=EXTENTS, ylim=EXTENTS)
        ax[1].set(xlim=[-1,1], ylim=[-1,1])

        # Do once: Place the colorbar on the second of the axes.
        if len(cbs) == 0:
            cbs.append(
                fig.colorbar(ax=ax[1], mappable=None, cmap=matplotlib.cm.afmhot, norm=NORM,
                             label=r'$\log_{10}$ M$_{\odot}$ kpc$^{-2}$')
            )

    # TODO: Parameterise this
    metadata = dict(title='NGC1904 Cartesian Projection')

    # Create a new writer. Use the file writer because the
    # stream writer seems to use too much memory.
    writer = FFMpegFileWriter(fps=24, metadata=metadata)
    with writer.saving(fig, 'movie.mp4', plt.rcParams['figure.dpi']):
        # Render each snapshot Frame as a movie frame.
        for frame in get_frames():
            render_frame(frame)
            writer.grab_frame()

        # Write the last frame out a few times.
        for _ in range(10):
            writer.grab_frame()
