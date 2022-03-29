import glob
import pandas as pd
import shutil
import os,sys

# moves avi files to local drive
def move(basepath,local_drive_path):
    """
    move all avi files from basepath (recursive) to local_drive_path
    and make csv of moved files containing their paths
    """
    # find avi files in all subfolders from basepath
    files = glob.glob(basepath + "/**/*.avi", recursive = True)
    # make and save csv of origin paths to the local drive
    df_ = pd.DataFrame()
    df_['files'] = files
    # section to not overwrite previous files
    if os.path.exists(os.path.join(local_drive_path,'files.csv')):
        df = pd.read_csv(os.path.join(local_drive_path,'files.csv'))
        df = pd.concat([df,df_],ignore_index=True)
        df = df.drop_duplicates(subset=['files'])
    else:
        df = df_    
    df.to_csv(os.path.join(local_drive_path,'files.csv'))
    # copy each video file to local drive
    for file in files:
        if not os.path.exists(os.path.join(local_drive_path,os.path.basename(file))):
            print('copy: ' + file)
            shutil.copyfile(file, os.path.join(local_drive_path,os.path.basename(file)))

# move csv files back to origin in basepath
def move_back(basepath,local_drive_path):
    """
    find deeplabcut results and copy them back to file paths in csv
    """
    # read files csv with origin basepaths
    df = pd.read_csv(os.path.join(local_drive_path,'files.csv'))
    # iter through each file
    for file in df.files:
        # check to see if csv exist for this video
        cur_file = os.path.splitext(os.path.basename(file))[0]
        cur_file_fullpath = glob.glob(os.path.join(local_drive_path,cur_file+"*.csv"), recursive = False)
        # if the csv exist, then move to origin basepath
        if len(cur_file_fullpath) > 0:
            destination = os.path.dirname(os.path.splitext(file)[0])
            if not os.path.exists(destination):
                print('WARNING: folder name change. can not copy results to '+destination)
                continue
            if not os.path.exists(os.path.join(destination,os.path.basename(cur_file_fullpath[0]))):

                print('copy: ' + cur_file_fullpath[0])
                shutil.copyfile(cur_file_fullpath[0], os.path.join(destination,os.path.basename(cur_file_fullpath[0])))
                
                # also move .h5 and .pickle
                h5_path = glob.glob(os.path.join(local_drive_path,cur_file+"*.h5"), recursive = False)
                shutil.copyfile(h5_path[0], os.path.join(destination,os.path.basename(h5_path[0])))

                pickle_path = glob.glob(os.path.join(local_drive_path,cur_file+"*.pickle"))
                shutil.copyfile(pickle_path[0], os.path.join(destination,os.path.basename(pickle_path[0])))

def main(basepath,local_drive_path):
    """
    run move and move_back
    """
    basepath = os.path.normpath(basepath)
    local_drive_path = os.path.normpath(local_drive_path)
    move(basepath,local_drive_path)
    move_back(basepath,local_drive_path)

# basepath = r'Z:\Data\Can\OLM21'
# local_drive_path = r'C:\Users\Cornell\dlc_videos_v2'

# basepath = os.path.normpath(basepath)
# local_drive_path = os.path.normpath(local_drive_path)
# main(basepath,local_drive_path)

# if __name__ == "__main__":
#     basepath = os.path.normpath(sys.argv[1])
#     local_drive_path = os.path.normpath(sys.argv[2])

#     print(basepath,local_drive_path)
#     main(basepath,local_drive_path)
