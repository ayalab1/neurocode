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

path_config_file = r"Y:\Praveen\social_plus_maze_multi_animal_dlc\social_plus_maze_multi_animal-ryanh-2022-04-23\config.yaml"

videos_folder = [r'Y:\Praveen\social_plus_maze_multi_animal_dlc\social_plus_maze_multi_animal-ryanh-2022-04-23\videos']

import sys
sys.path.append(r"C:\Users\praveen\Documents\GitHub\neurocode\behavior")
import find_and_move_videos
print('moving videos and results')
find_and_move_videos.main(r'Y:\Praveen\SocialBehavior',r'Y:\Praveen\social_plus_maze_multi_animal_dlc\social_plus_maze_multi_animal-ryanh-2022-04-23\videos')

deeplabcut.analyze_videos(path_config_file,videos_folder,shuffle=1,save_as_csv=True,dynamic=(True, .5, 10))

print('moving videos and results')
find_and_move_videos.main(r'Y:\Praveen\SocialBehavior',r'Y:\Praveen\social_plus_maze_multi_animal_dlc\social_plus_maze_multi_animal-ryanh-2022-04-23\videos')