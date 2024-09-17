#!/usr/bin/env python

"""
A script for turning NBODY6+P3T snapshots into a movie!

Expected input format is as follows (EBNF):

    FILE: SNAPSHOT*
    SNAPSHOT: HEADER "\n" (PARTICLE "\n")*
    HEADER: "T = " FLOAT
    PARTICLE: FLOAT{5}


This format matches the streaming output of the `stream` outputs program:

    $ stream <pos-file> | ./movie.py

The output of this script will be a 1280x720 MP4 @ 24fps.
"""

import collections
import sys
import numpy as np
import warnings
import matplotlib
import argparse
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegFileWriter
from matplotlib.path import Path

EX = Path([(-1, -1), (1, 1), (1, -1), (-1, 1)],
          [Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO])

NORM = matplotlib.colors.Normalize(0, 5)

# TODO: Read these from input
EXTENTS = [-10, 10]

class Frame:
    """
    A frame captures one snapshot from the input.
    """
    def __init__(self, t):
        """
        Initialise an empty frame with a specific timestamp.
        """
        self.t = t
        self.dat = None

    def __repr__(self):
        l = ' empty'
        if self.dat is not None:
            l = ''
        return f'<Frame {self.t:6.5}{l}>'

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
    r = np.array([ext, ext])
    m = frame.dat[:,0]

    # Convert positions to kpc
    if not gal:
        x = frame.dat[:,1] / 1000
        y = frame.dat[:,2] / 1000

    else:
        x = frame.dat[:,3] / 1000
        y = frame.dat[:,4] / 1000


    distance = np.abs(ext[1]-ext[0])
    h, _, _, i = ax.hist2d(x,
                           y,
                           bins=bins,
                           weights=m,
                           range=r,
                           cmap=matplotlib.cm.hot,
                           norm=NORM)
    sq=(distance/bins)**2
    d=h/sq

    a = np.ma.log10(d/1000).filled(0)
    i.set_array(a.T)
    return i

parser = argparse.ArgumentParser(
    prog='movie',
    description='Make a movie from your NBODY6+P3T output')
parser.add_argument('-t', '--title', type=str)
parser.add_argument('-s', '--start', help='The sim start time', type=float, default=9500.0)
parser.add_argument('-e', '--end', help='The sim end time', type=float, default=13500.0)

if __name__ == '__main__':
    args = parser.parse_args()

    px = 1/plt.rcParams['figure.dpi']

    plt.rcParams["figure.autolayout"] = True
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size'] = 18
    plt.style.use('dark_background')

    fig, ax = plt.subplots(ncols=2, figsize=(1920*px, 1080*px))
    #fig.tight_layout()

    ax[0].set(xlim=EXTENTS, ylim=EXTENTS)
    ax[0].set_xlabel('x [kpc]')
    ax[0].set_ylabel('y [kpc]')
    # Draw an X at the galaxy centre.
    cen = ax[0].plot([0], [0], 'kx', ms=15)

    # Write the time.
    txt = ax[0].text(-9.9, 10.5, ' ' * 15)

    ax[1].set_xlabel('x [kpc]')
    ax[1].set_ylabel('y [kpc]')
    meshes = []
    cbs = []

    def time(s):
        return float(s.rstrip().split()[-1]) - args.end

    def get_frames():
        """
        Stream the frames out one at a time.
        """
        assert(sys.stdin.readline().startswith('T'))
        frames = collections.deque([
            # Fix up the first frame time
            Frame(args.start-args.end)
        ])
        def snap():
            """
            Read each line from a frame until the next frame header
            is encountered.
            """
            for line in sys.stdin:
                if line.startswith('T'):
                    frames.append(Frame(time(line)))
                    return

                yield line

        # Read each frame from the input stream.
        while len(frames) > 0:
            frame = frames.popleft()
            frame.dat = np.loadtxt(snap())
            yield frame

    def render_frame(frame):
        print(frame)  # Report progress.
        txt.set_text(f'T = {frame.t:6.5} Myrs')

        # Remove old histograms from the plot.
        for mesh in meshes:
            mesh.remove()

        # Create new histograms.
        meshes[:] = [
            intensity(ax[0], frame, bins=1000, ext=EXTENTS),
            intensity(ax[1], frame, bins=250, ext=[-1, 1], gal=False)
        ]

        # Limits seem to need to be reset each time.
        ax[0].set(xlim=EXTENTS, ylim=EXTENTS)
        ax[1].set(xlim=[-1,1], ylim=[-1,1])

        # Do once: Place the colorbar on the second of the axes.
        if len(cbs) == 0:
            cbs.append(
                fig.colorbar(ax=ax[1],
                             mappable=None,
                             cmap=matplotlib.cm.afmhot,
                             norm=NORM,
                             label=r'Log M kpc$^{-2}$')
            )


    metadata = dict()
    if args.title:
        metadata['title'] = args.title

    # Create a new writer. Use the file writer because the
    # stream writer seems to use too much memory.
    writer = FFMpegFileWriter(fps=12.5, metadata=metadata)
    with writer.saving(fig, 'movie.mp4', plt.rcParams['figure.dpi']):
        # Render each snapshot Frame as a movie frame.
        for frame in get_frames():
            render_frame(frame)
            writer.grab_frame()

        # Write the last frame out a few times.
        for _ in range(15):
            writer.grab_frame()
