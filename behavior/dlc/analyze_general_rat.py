# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 16:33:27 2021

@author: Matt Isaacson
adapted by Ryan H
"""

import os

os.environ["DLClight"] = "False"
import tensorflow as tf

print(tf.__version__)

from tensorflow.python.client import device_lib

device_lib.list_local_devices()

os.environ["DLClight"] = "True"

import deeplabcut

path_config_file = r"D:\dlc_videos\general_rat_no_led\config.yaml"

# videos_folder = r'D:\dlc_videos\HMC2'
video_source = r"Z:\Data\HMC2"
import sys
import glob

# sys.path.append(r"D:\github\neurocode\behavior")
# import find_and_move_videos
print("moving videos and results")
# find_and_move_videos.main(video_source,videos_folder)
# find_and_move_videos.main(r'Y:\V1test\V1Jean',r'D:\dlc_videos\HMC2')

files = glob.glob(video_source + "/**/*.avi", recursive=True)

deeplabcut.analyze_videos(path_config_file, files, shuffle=1, save_as_csv=True)
deeplabcut.filterpredictions(path_config_file, files)

# print('moving videos and results')
# find_and_move_videos.main(video_source,videos_folder)
# find_and_move_videos.main(r'Y:\V1test\V1Jean',r'D:\dlc_videos\HMC2')
