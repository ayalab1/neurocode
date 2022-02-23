# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 16:25:22 2021

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
path_config_file = "D:/github/DeepLabCut/ojr-heath_ryan-2022-02-22/config.yaml"

videos_folder = ['C:/Users/Cornell/dlc_videos']

deeplabcut.analyze_videos(path_config_file,videos_folder,save_as_csv=True)