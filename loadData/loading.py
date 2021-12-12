import scipy.io as sio
import sys,os
import pandas as pd
import numpy as np
import glob
import warnings


def loadXML(path):
    """
    path should be the folder session containing the XML file
    Function returns :
        1. the number of channels
        2. the sampling frequency of the dat file or the eeg file depending of what is present in the folder
            eeg file first if both are present or both are absent
        3. the mappings shanks to channels as a dict
    Args:
        path : string
    Returns:
        int, int, dict

    by Guillaume Viejo    
    """
    if not os.path.exists(path):
        print("The path "+path+" doesn't exist; Exiting ...")
        sys.exit()
    listdir = os.listdir(path)
    xmlfiles = [f for f in listdir if f.endswith('.xml')]
    if not len(xmlfiles):
        print("Folder contains no xml files; Exiting ...")
        sys.exit()
    new_path = os.path.join(path, xmlfiles[0])

    from xml.dom import minidom
    xmldoc = minidom.parse(new_path)
    nChannels = xmldoc.getElementsByTagName('acquisitionSystem')[0].getElementsByTagName('nChannels')[0].firstChild.data
    fs_dat = xmldoc.getElementsByTagName('acquisitionSystem')[0].getElementsByTagName('samplingRate')[0].firstChild.data
    fs = xmldoc.getElementsByTagName('fieldPotentials')[0].getElementsByTagName('lfpSamplingRate')[0].firstChild.data

    shank_to_channel = {}
    groups = xmldoc.getElementsByTagName('anatomicalDescription')[0].getElementsByTagName('channelGroups')[0].getElementsByTagName('group')
    for i in range(len(groups)):
        shank_to_channel[i] = [int(child.firstChild.data) for child in groups[i].getElementsByTagName('channel')]
    return int(nChannels), int(fs), int(fs_dat), shank_to_channel

def loadLFP(path, n_channels=90, channel=64, frequency=1250.0, precision='int16'):
    if type(channel) is not list:
        f = open(path, 'rb')
        startoffile = f.seek(0, 0)
        endoffile = f.seek(0, 2)
        bytes_size = 2
        n_samples = int((endoffile-startoffile)/n_channels/bytes_size)
        duration = n_samples/frequency
        interval = 1/frequency
        f.close()
        with open(path, 'rb') as f:
            data = np.fromfile(f, np.int16).reshape((n_samples, n_channels))[:,channel]
            timestep = np.arange(0, len(data))/frequency
            # check if lfp time stamps exist
            lfp_ts_path = os.path.join(os.path.dirname(os.path.abspath(path)),'lfp_ts.npy')
            if os.path.exists(lfp_ts_path):
                timestep = np.load(lfp_ts_path).reshape(-1)

            return data, timestep # nts.Tsd(timestep, data, time_units = 's')
        
    elif type(channel) is list:
        f = open(path, 'rb')
        startoffile = f.seek(0, 0)
        endoffile = f.seek(0, 2)
        bytes_size = 2

        n_samples = int((endoffile-startoffile)/n_channels/bytes_size)
        duration = n_samples/frequency
        f.close()
        with open(path, 'rb') as f:
            data = np.fromfile(f, np.int16).reshape((n_samples, n_channels))[:,channel]
            timestep = np.arange(0, len(data))/frequency
            # check if lfp time stamps exist
            lfp_ts_path = os.path.join(os.path.dirname(os.path.abspath(path)),'lfp_ts.npy')
            if os.path.exists(lfp_ts_path):
                timestep = np.load(lfp_ts_path).reshape(-1)
            return data,timestep # nts.TsdFrame(timestep, data, time_units = 's')


def load_position(path,fs=39.0625):
    if not os.path.exists(path):
        print("The path "+path+" doesn't exist; Exiting ...")
        sys.exit()
    listdir = os.listdir(path)
    whlfiles = [f for f in listdir if f.endswith('.whl')]
    if not len(whlfiles):
        print("Folder contains no whl files; Exiting ...")
        sys.exit()
    new_path = os.path.join(path, whlfiles[0])
    df = pd.read_csv(new_path,delimiter="\t",header=0,names=['x1','y1','x2','y2'])
    df[df==-1] = np.nan
    return df,fs


def load_cell_metrics(basepath):
    """ 
    loader of cell-explorer cell_metrics.cellinfo.mat

    Inputs: basepath: path to folder with cell_metrics.cellinfo.mat
    outputs: df: data frame of single unit features
    data_: dict with data that does not fit nicely into a dataframe (waveforms, acgs, epochs, etc.)
    
    See https://cellexplorer.org/datastructure/standard-cell-metrics/ for details

    TODO: extract all fields from cell_metrics.cellinfo. There are more items that can be extracted

    - Ryan H
    """
    def extract_epochs(data):
        startTime = [ep['startTime'][0][0][0][0] for ep in data['cell_metrics']['general'][0][0]['epochs'][0][0][0]]
        stopTime = [ep['stopTime'][0][0][0][0] for ep in data['cell_metrics']['general'][0][0]['epochs'][0][0][0]]
        name = [ep['name'][0][0][0] for ep in data['cell_metrics']['general'][0][0]['epochs'][0][0][0]]

        epochs = pd.DataFrame()
        epochs['name'] = name
        epochs['startTime'] = startTime
        epochs['stopTime'] = stopTime
        return epochs

    def extract_general(data):
        # extract fr per unit with lag zero to ripple
        try:
            ripple_fr = [ev.T[0] for ev in data['cell_metrics']['events'][0][0]['ripples'][0][0][0]]
        except:
            ripple_fr = []
        # extract spikes times
        spikes = [spk.T[0] for spk in data['cell_metrics']['spikes'][0][0]['times'][0][0][0]]
        # extract epochs
        try:
            epochs = extract_epochs(data)
        except:
            epochs = []
        # extract avg waveforms (one wavefrom per channel on shank)
        try:
            waveforms = [w.T for w in data['cell_metrics']['waveforms'][0][0][0][0][0][0]]
        except:
            waveforms = [w.T for w in data['cell_metrics']['waveforms'][0][0][0]]
        # extract chanCoords
        try:
            chanCoords_x = data['cell_metrics']['general'][0][0]['chanCoords'][0][0][0][0]['x'].T[0]
            chanCoords_y = data['cell_metrics']['general'][0][0]['chanCoords'][0][0][0][0]['y'].T[0]
        except:
            chanCoords_x = []
            chanCoords_y = []

        # add to dictionary 
        data_ = {
            "acg_wide": data['cell_metrics']['acg'][0][0]['wide'][0][0],
            "acg_narrow": data['cell_metrics']['acg'][0][0]['narrow'][0][0],
            "acg_log10": data['cell_metrics']['acg'][0][0]['log10'][0][0],
            "ripple_fr": ripple_fr,
            "chanCoords_x": chanCoords_x,
            "chanCoords_y": chanCoords_y,
            "epochs": epochs,
            "spikes": spikes,
            "waveforms": waveforms
            }
        return data_ 

    def un_nest_df(df):
        # Un-nest some strings are nested within brackets (a better solution exists...)
        # locate and iterate objects in df
        for item in df.keys()[df.dtypes =="object"]:
            # if you can get the size of the first item with [0], it is nested
            # otherwise it fails and is not nested
            try:
                df[item][0][0].size
                # the below line is from: https://www.py4u.net/discuss/140913
                df[item] = df[item].str.get(0)
            except:
                continue
        return df

    filename = glob.glob(os.path.join(basepath,'*.cell_metrics.cellinfo.mat'))[0]

    # check if saved file exists
    if not os.path.exists(filename):
        warnings.warn("file does not exist")
        return 

    # load cell_metrics file
    data = sio.loadmat(filename)

    # construct data frame with features per neuron
    df = pd.DataFrame()
    # count units
    n_cells = data['cell_metrics']['UID'][0][0][0].size
    dt = data['cell_metrics'].dtype
    for dn in dt.names:
        # check if var has the right n of units and is a vector
        try:
            if ((data['cell_metrics'][dn][0][0][0][0].size == 1) &
                    (data['cell_metrics'][dn][0][0][0].size == n_cells)):
                    
                df[dn] = data['cell_metrics'][dn][0][0][0]
        except:
            continue
        
    # add column for bad label tag    
    try:
        bad_units = data['cell_metrics']['tags'][0][0]['Bad'][0][0][0]
        df['bad_unit'] = [False]*df.shape[0]
        for uid in bad_units:
            df.loc[df.UID == uid,'bad_unit'] = True
    except:
        df['bad_unit'] = [False]*df.shape[0]  

    # add data from general metrics        
    df['basename'] = data['cell_metrics']['general'][0][0]['basename'][0][0][0]
    df['basepath'] = data['cell_metrics']['general'][0][0]['basepath'][0][0][0]
    df['sex'] = data['cell_metrics']['general'][0][0]['animal'][0][0]['sex'][0][0][0]
    df['species'] = data['cell_metrics']['general'][0][0]['animal'][0][0]['species'][0][0][0]
    df['strain'] = data['cell_metrics']['general'][0][0]['animal'][0][0]['strain'][0][0][0]
    df['geneticLine'] = data['cell_metrics']['general'][0][0]['animal'][0][0]['geneticLine'][0][0][0]
    df['cellCount'] = data['cell_metrics']['general'][0][0]['cellCount'][0][0][0][0]

    # fix nesting issue for strings
    df = un_nest_df(df)

    # extract other general data and put into dict    
    data_ = extract_general(data)

    return df,data_

def load_SWRunitMetrics(basepath):
    """
    load_SWRunitMetrics loads SWRunitMetrics.mat into pandas dataframe

    returns pandas dataframe with the following fields
        particip: the probability of participation into ripples for each unit
        FRall: mean firing rate during ripples
        FRparticip: mean firing rate for ripples with at least 1 spike
        nSpkAll: mean number of spikes in all ripples
        nSpkParticip: mean number of spikes in ripples with at least 1 spike
        epoch: behavioral epoch label
    """

    def extract_swr_epoch_data(data,epoch):
        # get var names
        dt = data['SWRunitMetrics'][epoch][0][0].dtype

        df2 = pd.DataFrame()

        # get n units
        # there might be other fields within here like the epoch timestamps
        # skip those by returning empty df
        try:
            n_cells = data['SWRunitMetrics'][epoch][0][0][0]['particip'][0].shape[0]
        except:
            return df2

        for dn in dt.names:
            if (
                (data['SWRunitMetrics'][epoch][0][0][0][dn][0].shape[1] == 1) &
                (data['SWRunitMetrics'][epoch][0][0][0][dn][0].shape[0] == n_cells)
                ):
                df2[dn] = data['SWRunitMetrics'][epoch][0][0][0][dn][0].T[0]
        df2['epoch'] = epoch
        return df2

    filename = glob.glob(os.path.join(basepath,'*.SWRunitMetrics.mat'))[0]

    # check if saved file exists
    if not os.path.exists(filename):
        warnings.warn("file does not exist")
        return pd.DataFrame()

    # load file
    data = sio.loadmat(filename)

    df2 = pd.DataFrame()
    # loop through each available epoch and pull out contents
    for epoch in data['SWRunitMetrics'].dtype.names:
        if data['SWRunitMetrics'][epoch][0][0].size>0: # not empty
            
            # call content extractor 
            df_ = extract_swr_epoch_data(data,epoch)

            # append conents to overall data frame
            if df_.size>0:
                df2 = df2.append(df_,ignore_index=True)

    return df2

def load_ripples_events(basepath):
    """
    load info from ripples.events.mat and store within df

    basepath: path to your session where ripples.events.mat is
    
    returns pandas dataframe with the following fields
        start: start time of ripple
        stop: end time of ripple
        peaks: peak time of ripple
        amplitude: envlope value at peak time
        duration: ripple duration
        frequency: insta frequency at peak
        detectorName: the name of ripple detector used
        event_spk_thres: 1 or 0 for if a mua thres was used 
        basepath: path name
        basename: session id
        animal: animal id *

        * Note that basepath/basename/animal relies on specific folder 
        structure and may be incorrect for some data structures
    """

    # locate .mat file
    filename = glob.glob(basepath+os.sep+'*ripples.events.mat')[0]

    # check if saved file exists
    if not os.path.exists(filename):
        warnings.warn("file does not exist")
        return pd.DataFrame()

    # load matfile
    data = sio.loadmat(filename)

    # make data frame of known fields 
    df = pd.DataFrame()
    try:
        df['start'] = data['ripples']['timestamps'][0][0][:,0]
        df['stop'] = data['ripples']['timestamps'][0][0][:,1]
    except:
        df['start'] = data['ripples']['times'][0][0][:,0]
        df['stop'] = data['ripples']['times'][0][0][:,1]

    df['peaks'] = data['ripples']['peaks'][0][0]
    try:
        df['amplitude'] = data['ripples']['amplitude'][0][0]
        df['duration'] = data['ripples']['duration'][0][0]
        df['frequency'] = data['ripples']['frequency'][0][0]
    except:
        df['amplitude'] = np.nan
        df['duration'] = np.nan
        df['frequency'] = np.nan
     
    try:
        df['detectorName'] = data['ripples']['detectorinfo'][0][0]['detectorname'][0][0][0]
    except:
        df['detectorName'] = data['ripples']['detectorName'][0][0][0]

    # find ripple channel (this can be in several places depending on the file)
    try:
        df['ripple_channel'] = data['ripples']['detectorinfo'][0][0]['detectionparms'][0][0]['Channels'][0][0][0][0]
    except:
        try:
            df['ripple_channel'] = data['ripples']['detectorParams'][0][0]['channel'][0][0][0][0]
        except:
            try:
                df['ripple_channel'] = data['ripples']['detectorinfo'][0][0]['detectionparms'][0][0]['channel'][0][0][0][0]
            except:
                df['ripple_channel'] = data['ripples']['detectorinfo'][0][0]['detectionparms'][0][0]['ripple_channel'][0][0][0][0]


    dt = data['ripples'].dtype
    if "eventSpikingParameters" in dt.names:
        df['event_spk_thres'] = 1
    else:
        df['event_spk_thres'] = 0

    # get basename and animal
    normalized_path = os.path.normpath(filename)
    path_components = normalized_path.split(os.sep)
    df['basepath'] = basepath  
    df['basename'] = path_components[-2]
    df['animal'] = path_components[-3]

    return df


def load_theta_rem_shift(basepath):
    """
    load_theta_rem_shift: loads matlab structure from get_rem_shift.m
    """
    try:
        filename = glob.glob(basepath+os.sep+'*theta_rem_shift.mat')[0]
    except:
        # warnings.warn("file does not exist")
        return pd.DataFrame(),np.nan

    # check if saved file exists
    if not os.path.exists(filename):
        # warnings.warn("file does not exist")
        return pd.DataFrame(),np.nan

    data = sio.loadmat(filename)

    df = pd.DataFrame()

    df["circ_dist"] = data['rem_shift_data']['circ_dist'][0][0][0]
    df["rem_shift"] = data['rem_shift_data']['rem_shift'][0][0][0]
    df["non_rem_shift"] = data['rem_shift_data']['non_rem_shift'][0][0][0]
    
    # rem metrics
    df["m_rem"] = data['rem_shift_data']['PhaseLockingData_rem'][0][0]['phasestats'][0][0]['m'][0][0][0]
    df["r_rem"] = data['rem_shift_data']['PhaseLockingData_rem'][0][0]['phasestats'][0][0]['r'][0][0][0]
    df["k_rem"] = data['rem_shift_data']['PhaseLockingData_rem'][0][0]['phasestats'][0][0]['k'][0][0][0]
    df["p_rem"] = data['rem_shift_data']['PhaseLockingData_rem'][0][0]['phasestats'][0][0]['p'][0][0][0]
    df["mode_rem"] = data['rem_shift_data']['PhaseLockingData_rem'][0][0]['phasestats'][0][0]['mode'][0][0][0]
    
    # wake metrics
    df["m_wake"] = data['rem_shift_data']['PhaseLockingData_wake'][0][0]['phasestats'][0][0]['m'][0][0][0]
    df["r_wake"] = data['rem_shift_data']['PhaseLockingData_wake'][0][0]['phasestats'][0][0]['r'][0][0][0]
    df["k_wake"] = data['rem_shift_data']['PhaseLockingData_wake'][0][0]['phasestats'][0][0]['k'][0][0][0]
    df["p_wake"] = data['rem_shift_data']['PhaseLockingData_wake'][0][0]['phasestats'][0][0]['p'][0][0][0]
    df["mode_wake"] = data['rem_shift_data']['PhaseLockingData_wake'][0][0]['phasestats'][0][0]['mode'][0][0][0]


    def get_distros(data,state):
        return np.vstack(data['rem_shift_data'][state][0][0]['phasedistros'][0][0].T)

    def get_spikephases(data,state):
        return data['rem_shift_data'][state][0][0]['spkphases'][0][0][0]

    # add to dictionary 
    data_dict = {
                    "rem": 
                    {
                        "phasedistros": get_distros(data,'PhaseLockingData_rem'),
                        "spkphases":get_spikephases(data,'PhaseLockingData_rem')
                    },
                    "wake":
                    {
                        'phasedistros': get_distros(data,'PhaseLockingData_wake'),
                        "spkphases":get_spikephases(data,'PhaseLockingData_wake')
                    }
                }

    return df,data_dict


def load_SleepState_states(basepath):
    """ 
    loader of SleepState.states.mat

    returns dict of structures contents.

    TODO: extract more from file, this extracts the basics for now.

    """
    filename = glob.glob(os.path.join(basepath,'*.SleepState.states.mat'))[0]
    
    # check if saved file exists
    if not os.path.exists(filename):
        warnings.warn("file does not exist")
        return 

    # load cell_metrics file
    data = sio.loadmat(filename)

    # get epoch id
    wake_id = np.where(data["SleepState"]["idx"][0][0]["statenames"][0][0][0] == "WAKE")[0][0]+1
    rem_id = np.where(data["SleepState"]["idx"][0][0]["statenames"][0][0][0] == "REM")[0][0]+1
    nrem_id = np.where(data["SleepState"]["idx"][0][0]["statenames"][0][0][0] == "NREM")[0][0]+1

    # get states and timestamps vectors
    states = data["SleepState"]["idx"][0][0]["states"][0][0]
    timestamps = data["SleepState"]["idx"][0][0]["timestamps"][0][0]

    # set up dict
    dict_ = {"wake_id": wake_id, "rem_id": rem_id, "nrem_id": nrem_id,
        "states": states, "timestamps": timestamps}

    # iter through states and add to dict   
    dt = data["SleepState"]["ints"][0][0].dtype
    for dn in dt.names:
        dict_[dn] = data["SleepState"]["ints"][0][0][dn][0][0]

    return dict_

def load_animal_behavior(basepath):
    """
    load_animal_behavior loads basename.animal.behavior.mat files created by general_behavior_file.m
    The output is a pandas data frame with [time,x,y,z,linerized,speed,acceleration,trials,epochs]

    Ryan H 2021
    """

    filename = glob.glob(os.path.join(basepath,'*.animal.behavior.mat'))[0]

    # check if saved file exists
    if not os.path.exists(filename):
        warnings.warn("file does not exist")
        return 

    def extract_epochs(data):
        startTime = [ep['startTime'][0][0][0][0] for ep in data['behavior']['epochs'][0][0][0] if len(ep[0]) > 0]
        stopTime = [ep['stopTime'][0][0][0][0] for ep in data['behavior']['epochs'][0][0][0] if len(ep[0]) > 0]
        name = [ep['name'][0][0][0] for ep in data['behavior']['epochs'][0][0][0] if len(ep[0]) > 0]

        epochs = pd.DataFrame()
        epochs['name'] = name
        epochs['startTime'] = startTime
        epochs['stopTime'] = stopTime
        return epochs

    # load cell_metrics file
    data = sio.loadmat(filename)

    trials = data['behavior']['trials'][0][0]
    
    epochs = extract_epochs(data)

    df = pd.DataFrame()
    try:
        df['time'] = data['behavior']['time'][0][0][0]
    except:
        warnings.warn("no tracking data")
        return pd.DataFrame()
        
    try:
        df['x'] = data['behavior']['position'][0][0]['x'][0][0][0]
    except:
        df['x'] = np.nan
    try:
        df['y'] = data['behavior']['position'][0][0]['y'][0][0][0]
    except:
        df['y'] = np.nan
    try:
        df['z'] = data['behavior']['position'][0][0]['z'][0][0][0]
    except:
        df['z'] = np.nan
    try:
        df['linerized'] = data['behavior']['position'][0][0]['linerized'][0][0][0]
    except:
        df['linerized'] = np.nan

    df['speed'] = data['behavior']['speed'][0][0][0]
    df['acceleration'] = data['behavior']['acceleration'][0][0][0]

    for t in range(trials.shape[0]):
        idx = (df.time >= trials[t,0]) & (df.time <= trials[t,1])
        df.loc[idx,'trials'] = t

    for t in range(epochs.shape[0]):
        idx = (df.time >= epochs.startTime.iloc[t]) & (df.time <= epochs.stopTime.iloc[t])
        df.loc[idx,'epochs'] = epochs.name.iloc[t] 

    return df

def load_epoch(basepath):
    """
    Loads epoch info from cell explorer basename.session and stores in df
    """

    filename = glob.glob(os.path.join(basepath,'*.session.mat'))[0]

    # check if saved file exists
    if not os.path.exists(filename):
        warnings.warn("file does not exist")
        return pd.DataFrame()

    # load file
    data = sio.loadmat(filename)

    name = []
    df_temp = pd.DataFrame()
    df_save = pd.DataFrame()
    for epoch in data['session']['epochs'][0][0][0]:
        if len(epoch[0]) == 0:
            continue
        dt = epoch[0].dtype
        for dn in dt.names:
            try:
                df_temp[dn] = epoch[0][0][dn][0]
            except:
                df_temp[dn] = ""
        df_save = df_save.append(df_temp,ignore_index=True)
        name.append(epoch[0]['name'][0][0])
    df_save['name'] = name

    return df_save