# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 20:40:34 2019

@author: hindesa
"""

import plotly.graph_objs as go
import plotly.offline as py

import numpy as np
from ipywidgets import interactive, HBox, VBox

x = y = np.arange(-5, 5, 0.1)
yt = x[:, np.newaxis]
z = np.cos(x * yt) + np.sin(x * yt) * 2

f = go.FigureWidget(
    data=[
        go.Surface(z=z, x=x, y=y,
                   colorscale='Viridis')],
    layout=go.Layout(scene=go.layout.Scene(
        camera=go.layout.scene.Camera(
            up=dict(x=0, y=0, z=1),
            center=dict(x=0, y=0, z=0),
            eye=dict(x=1.25, y=1.25, z=1.25))
    ))
)


def update_z(frequency):
    f.data[0].z = np.cos(x * yt * frequency / 10.0) + np.sin(x * yt * frequency / 10.0) * 2


freq_slider = interactive(update_z, frequency=(1, 50, 0.1))
vb = VBox((f, freq_slider))
vb.layout.align_items = 'center'
vb