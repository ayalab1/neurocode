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

from tensorflow.pyzthon.client import device_lib
device_lib.list_local_devices()

os.environ["DLClight"] = "True"

import deeplabcut

path_config_file = r"D:\behavior\Peacemaze-Hongyu-2022-06-16\config.yaml"

videos_folder = [r'D:\behavior\dlc_video']

import sys
sys.path.append(r"D:\github\neurocode\behavior")
import find_and_move_videos
print('moving videos and results')
find_and_move_videos.main(r'N:\Hongyu\maze1ormaze2\HC1',r'D:\behavior\dlc_video')
#find_and_move_videos.main(r'Y:\V1test\V1Jean',r'D:\dlc_videos\v1_implanted')

deeplabcut.analyze_videos(path_config_file,videos_folder,shuffle=1,save_as_csv=True,dynamic=(True, .5, 10))

print('moving videos and results')
find_and_move_videos.main(r'N:\Hongyu\maze1ormaze2\HC1',r'D:\behavior\dlc_video')
#find_and_move_videos.main(r'Y:\V1test\V1Jean',r'D:\dlc_videos\v1_implanted')
