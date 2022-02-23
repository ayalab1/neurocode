# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 16:33:27 2021

@author: Matt Isaacson
"""

import os
os.environ["DLClight"] = "False"
import tensorflow as tf
print(tf.__version__)

from tensorflow.python.client import device_lib
device_lib.list_local_devices()

os.environ["DLClight"] = "True"

import deeplabcut
path_config_file = r"C:/Users/Schaf/Documents/DLC Networks/YMaze-Matt-2021-11-22/config.yaml"

videos_folder = ['C:/Users/Schaf/Desktop/Y-M Videos to Analyze']

deeplabcut.analyze_videos(path_config_file,videos_folder)