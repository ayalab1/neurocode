from logging import exception
import scipy.io as sio
import sys,os
import pandas as pd
import numpy as np
import glob
import nelpy as nel
import warnings
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

def loadXML(basepath):
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
    # check if saved file exists
    try:
        basename = os.path.basename(basepath)
        filename = glob.glob(os.path.join(basepath,basename+'.xml'))[0]
    except:
        warnings.warn("file does not exist")
        return 

    from xml.dom import minidom
    xmldoc = minidom.parse(filename)
    nChannels = xmldoc.getElementsByTagName('acquisitionSystem')[0].getElementsByTagName('nChannels')[0].firstChild.data
    fs_dat = xmldoc.getElementsByTagName('acquisitionSystem')[0].getElementsByTagName('samplingRate')[0].firstChild.data
    fs = xmldoc.getElementsByTagName('fieldPotentials')[0].getElementsByTagName('lfpSamplingRate')[0].firstChild.data

    shank_to_channel = {}
    groups = (
        xmldoc.getElementsByTagName("anatomicalDescription")[0]
        .getElementsByTagName("channelGroups")[0]
        .getElementsByTagName("group")
    )
    for i in range(len(groups)):
        shank_to_channel[i] = [int(child.firstChild.data) for child in groups[i].getElementsByTagName('channel')]
    return int(nChannels), int(fs), int(fs_dat), shank_to_channel

def loadLFP(basepath, n_channels=90, channel=64, frequency=1250.0, precision='int16',ext='lfp'):
    if ext == 'lfp':
        try:
            path = glob.glob(os.path.join(basepath,os.path.basename(basepath)+'*.lfp'))[0]
        except:
            path = glob.glob(os.path.join(basepath,os.path.basename(basepath)+'*.eeg'))[0]
    if ext == 'dat':
            path = glob.glob(os.path.join(basepath,os.path.basename(basepath)+'*.dat'))[0]

    # check if saved file exists
    if not os.path.exists(path):
        warnings.warn("file does not exist")
        return 

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

def load_position(basepath,fs=39.0625):
    if not os.path.exists(basepath):
        print("The path "+basepath+" doesn't exist; Exiting ...")
        sys.exit()
    listdir = os.listdir(basepath)
    whlfiles = [f for f in listdir if f.endswith('.whl')]
    if not len(whlfiles):
        print("Folder contains no whl files; Exiting ...")
        sys.exit()
    new_path = os.path.join(basepath, whlfiles[0])
    df = pd.read_csv(new_path,delimiter="\t",header=0,names=['x1','y1','x2','y2'])
    df[df==-1] = np.nan
    return df,fs

def load_all_cell_metrics(basepaths):
    """
    load cell metrics from multiple sessions

    Input: 
            basepaths: list of basepths, can be pandas column
    Output:
            cell_metrics: concatenated pandas df with metrics

    Note: to get waveforms, spike times, etc. use load_cell_metrics
    """
    import multiprocessing
    from joblib import Parallel, delayed

    # to speed up, use parallel
    num_cores = multiprocessing.cpu_count()   
    cell_metrics = Parallel(n_jobs=num_cores)(
            delayed(load_cell_metrics)(basepath,True) for basepath in basepaths
        )

    return pd.concat(cell_metrics,ignore_index=True)


def load_cell_metrics(basepath,only_metrics=False):
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
        # extract avg waveforms
        try:
            waveforms = np.vstack(data['cell_metrics']['waveforms'][0][0]["filt"][0][0][0])
        except:
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
    # have dedicated var as this tag is important   
    try:
        bad_units = data['cell_metrics']['tags'][0][0]['Bad'][0][0][0]
        df['bad_unit'] = [False]*df.shape[0]
        for uid in bad_units:
            df.loc[df.UID == uid,'bad_unit'] = True
    except:
        df['bad_unit'] = [False]*df.shape[0]  

    # load in tag
    try:
        dt = data['cell_metrics']['tags'][0][0].dtype
        if len(dt) > 0:
            # iter through each tag
            for dn in dt.names:
                # set up column for tag
                df['tags_'+dn] = [False]*df.shape[0]
                # iter through uid 
                for uid in data['cell_metrics']['tags'][0][0][dn][0][0][0]:
                    df.loc[df.UID == uid,'tags_'+dn] = True 
    except:
        pass
    
    # add data from general metrics        
    df['basename'] = data['cell_metrics']['general'][0][0]['basename'][0][0][0]
    df['basepath'] = basepath
    df['sex'] = data['cell_metrics']['general'][0][0]['animal'][0][0]['sex'][0][0][0]
    df['species'] = data['cell_metrics']['general'][0][0]['animal'][0][0]['species'][0][0][0]
    df['strain'] = data['cell_metrics']['general'][0][0]['animal'][0][0]['strain'][0][0][0]
    try:
        df['geneticLine'] = data['cell_metrics']['general'][0][0]['animal'][0][0]['geneticLine'][0][0][0]
    except:
        pass
    df['cellCount'] = data['cell_metrics']['general'][0][0]['cellCount'][0][0][0][0]

    # fix nesting issue for strings
    df = un_nest_df(df)

    # convert nans within tags columns to false
    cols = df.filter(regex='tags_').columns
    df[cols] = df[cols].replace({np.nan:False})
    
    if only_metrics:
        return df

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

    try:
        filename = glob.glob(os.path.join(basepath,'*.SWRunitMetrics.mat'))[0]
    except:
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
                df2 = pd.concat([df2,df_],ignore_index=True)

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
    try:
        filename = glob.glob(basepath+os.sep+'*ripples.events.mat')[0]
    except:
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
    except:
        df['amplitude'] = np.nan
    try:
        df['duration'] = data['ripples']['duration'][0][0]
    except:
        df['duration'] = np.nan
    try:
        df['frequency'] = data['ripples']['frequency'][0][0]
    except:
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

def load_dentate_spike(basepath):
    """
    load info from DS*.events.mat and store within df
    basepath: path to your session where DS*.events.mat is

    returns pandas dataframe with the following fields
        start: start time of DS
        stop: end time of DS
        peaks: peak time of DS
        amplitude: envlope value at peak time
        duration: DS duration
        detectorName: the name of DS detector used
        basepath: path name
        basename: session id
        animal: animal id *
        * Note that basepath/basename/animal relies on specific folder
        structure and may be incorrect for some data structures
    """

    def extract_data(s_type, data):
        # make data frame of known fields
        df = pd.DataFrame()
        df["start"] = data[s_type]["timestamps"][0][0][:, 0]
        df["stop"] = data[s_type]["timestamps"][0][0][:, 1]
        df["peaks"] = data[s_type]["peaks"][0][0]
        df["event_label"] = s_type
        df["amplitude"] = data[s_type]["amplitudes"][0][0]
        df["duration"] = data[s_type]["duration"][0][0]
        df["amplitudeUnits"] = data[s_type]["amplitudeUnits"][0][0][0]
        df["detectorName"] = data[s_type]["detectorinfo"][0][0]["detectorname"][0][0][0]
        df["ml_channel"] = data[s_type]["detectorinfo"][0][0]["ml_channel"][0][0][0][0]
        df["h_channel"] = data[s_type]["detectorinfo"][0][0]["h_channel"][0][0][0][0]
        return df

    # locate .mat file
    df = pd.DataFrame()
    for s_type in ["DS1", "DS2"]:
        filename = glob.glob(basepath + os.sep + "*" + s_type + ".events.mat")[0]
        # load matfile
        data = sio.loadmat(filename)
        # pull out data
        df = pd.concat([df,extract_data(s_type, data)], ignore_index=True)

    # get basename and animal
    normalized_path = os.path.normpath(filename)
    path_components = normalized_path.split(os.sep)
    df["basepath"] = basepath
    df["basename"] = path_components[-2]
    df["animal"] = path_components[-3]

    return df
    
def load_theta_rem_shift(basepath):
    """
    load_theta_rem_shift: loads matlab structure from get_rem_shift.m
    """
    try:
        filename = glob.glob(basepath+os.sep+'*theta_rem_shift.mat')[0]
    except:
        warnings.warn("file does not exist")
        return pd.DataFrame(),np.nan

    data = sio.loadmat(filename)

    df = pd.DataFrame()

    df["UID"] = data['rem_shift_data']['UID'][0][0][0]
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
    try:
        filename = glob.glob(os.path.join(basepath,'*.SleepState.states.mat'))[0]
    except:
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
    try:
        filename = glob.glob(os.path.join(basepath,'*.animal.behavior.mat'))[0]
    except:
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
        df['time'] = data['behavior']['timestamps'][0][0][0]
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
        df['linearized'] = data['behavior']['position'][0][0]['linearized'][0][0][0]
    except:
        df['linearized'] = np.nan
    try:
        df['speed'] = data['behavior']['speed'][0][0][0]
    except:
        df['speed'] = np.nan
    try:
        df['acceleration'] = data['behavior']['acceleration'][0][0][0]
    except:
        df['acceleration'] = np.nan

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
    try:
        filename = glob.glob(os.path.join(basepath,'*.session.mat'))[0]
    except:
        warnings.warn("file does not exist")
        return pd.DataFrame()

    # load file
    data = sio.loadmat(filename)

    name = []
    environment = []
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
        df_save = pd.concat([df_save,df_temp],ignore_index=True)
        name.append(epoch[0]['name'][0][0])
        try:
            environment.append(epoch[0]['environment'][0][0])
        except:
            environment.append('unknown')
    df_save['name'] = name
    df_save['environment'] = environment

    return df_save
    
def load_brain_regions(basepath):
    """
    Loads brain region info from cell explorer basename.session and stores in dict

    Example:
        Input:
            brainRegions = load_epoch("Z:\\Data\\GirardeauG\\Rat09\\Rat09-20140327")
            print(brainRegions.keys())
            print(brainRegions['CA1'].keys())
            print(brainRegions['CA1']['channels'])
            print(brainRegions['CA1']['electrodeGroups'])
        output:
            dict_keys(['CA1', 'Unknown', 'blv', 'bmp', 'ven'])
            dict_keys(['channels', 'electrodeGroups'])
            [145 146 147 148 149 153 155 157 150 151 154 159 156 152 158 160 137 140
            129 136 138 134 130 132 142 143 144 141 131 139 133 135]
            [17 18 19 20]
    """
    filename = glob.glob(os.path.join(basepath, "*.session.mat"))[0]

    # load file
    data = sio.loadmat(filename)
    data = data["session"]
    brainRegions = {}
    for dn in data["brainRegions"][0][0].dtype.names:
        brainRegions[dn] = {
            "channels": data["brainRegions"][0][0][dn][0][0][0][0][0][0],
            "electrodeGroups": data["brainRegions"][0][0][dn][0][0][0][0][1][0],
        }

    return brainRegions

def get_animal_id(basepath):
    """ return animal ID from basepath using basename.session.mat"""
    try:
        filename = glob.glob(os.path.join(basepath,'*.session.mat'))[0]
    except:
        warnings.warn("file does not exist")
        return pd.DataFrame()

    # load file
    data = sio.loadmat(filename)
    return data['session'][0][0]['animal'][0][0]['name'][0]

def load_basic_data(basepath):
    nChannels, fs, fs_dat, shank_to_channel = loadXML(basepath)
    ripples = load_ripples_events(basepath)
    cell_metrics,data = load_cell_metrics(basepath)
    return cell_metrics,data,ripples,fs_dat

def load_spikes(basepath,
                putativeCellType=[], # restrict spikes to putativeCellType
                brainRegion=[], # restrict spikes to brainRegion
                bad_unit=False, # false for not loading bad cells
                brain_state=[], # restrict spikes to brainstate
                other_metric=None, # restrict spikes to other_metric
                other_metric_value=None # restrict spikes to other_metric_value
                ):
    """ 
    Load specific cells' spike times
    """
    if not isinstance(putativeCellType, list):
        putativeCellType = [putativeCellType]
    if not isinstance(brainRegion, list):
        brainRegion = [brainRegion]

    _,_,fs_dat,_ = loadXML(basepath)

    cell_metrics,data = load_cell_metrics(basepath)

    st = np.array(data['spikes'],dtype=object)

    # restrict cell metrics                      
    if len(putativeCellType) > 0:
        restrict_idx = []
        for cell_type in putativeCellType:
            restrict_idx.append(cell_metrics.putativeCellType.str.contains(cell_type).values) 
        restrict_idx = np.any(restrict_idx,axis=0)
        cell_metrics = cell_metrics[restrict_idx]
        st = st[restrict_idx]

    if len(brainRegion) > 0:
        restrict_idx = []
        for brain_region in brainRegion:
            restrict_idx.append(cell_metrics.brainRegion.str.contains(brain_region).values) 
        restrict_idx = np.any(restrict_idx,axis=0)
        cell_metrics = cell_metrics[restrict_idx]
        st = st[restrict_idx]

    # restrict cell metrics by arbitrary metric
    if other_metric is not None:
        # make other_metric_value a list if not already
        if not isinstance(other_metric, list):
            other_metric = [other_metric]
        if not isinstance(other_metric_value, list):
            other_metric_value = [other_metric_value]
        # check that other_metric_value is the same length as other_metric
        if len(other_metric) != len(other_metric_value):
            raise ValueError('other_metric and other_metric_value must be of same length')

        restrict_idx = []
        for metric,value in zip(other_metric,other_metric_value):
            restrict_idx.append(cell_metrics[metric].str.contains(value).values) 
        restrict_idx = np.any(restrict_idx,axis=0)
        cell_metrics = cell_metrics[restrict_idx]
        st = st[restrict_idx]

    restrict_idx = cell_metrics.bad_unit.values==bad_unit
    cell_metrics = cell_metrics[restrict_idx]
    st = st[restrict_idx]

    # get spike train array
    try:
        st = nel.SpikeTrainArray(timestamps=st, fs=fs_dat)
    except: # if only single cell... should prob just skip session
        st = nel.SpikeTrainArray(timestamps=st[0], fs=fs_dat)

    if len(brain_state) > 0:
        # get brain states        
        brain_states = ['WAKEstate', 'NREMstate', 'REMstate', 'THETA', 'nonTHETA']
        if brain_state not in brain_states:
            assert print('not correct brain state. Pick one',brain_states) 
        else:                                    
            state_dict = load_SleepState_states(basepath)
            state_epoch = nel.EpochArray(state_dict[brain_state])
            st = st[state_epoch]

    return st,cell_metrics

def load_deepSuperficialfromRipple(basepath,bypass_mismatch_exception=False):
    """
    Load deepSuperficialfromRipple file created by classification_DeepSuperficial.m
    
    """
    # locate .mat file
    file_type = "*.deepSuperficialfromRipple.channelinfo.mat"
    filename = glob.glob(basepath + os.sep + file_type)[0]

    # load matfile
    data = sio.loadmat(filename)

    channel_df = pd.DataFrame()
    name = "deepSuperficialfromRipple"

    # sometimes more channels positons will be in deepSuperficialfromRipple than in xml
    #   this is because they used channel id as an index. 
    channel_df = pd.DataFrame()
    channels = np.hstack(data[name]["channel"][0][0])*np.nan
    shanks = np.hstack(data[name]["channel"][0][0])*np.nan

    channels_, shanks_ = zip(
        *[(values[0], np.tile(shank, len(values[0]))) for shank, values in enumerate(data[name]["ripple_channels"][0][0][0])]
    )
    channel_sort_idx = np.hstack(channels_)-1
    channels[channel_sort_idx] = np.hstack(channels_)
    shanks[channel_sort_idx] = np.hstack(shanks_) + 1

    channel_df["channel"] = channels
    channel_df.loc[np.arange(len(channel_sort_idx)),"channel_sort_idx"] = channel_sort_idx
    channel_df["shank"] = shanks

    # add distance from pyr layer (will only be accurate if polarity rev)
    channel_df["channelDistance"] = data[name]["channelDistance"][0][0].T[0]

    # add channel class (deep or superficial)
    channelClass = []
    for item in data[name]["channelClass"][0][0]:
        try:
            channelClass.append(item[0][0])
        except:
            channelClass.append("unknown")
    channel_df["channelClass"] = channelClass

    # add if shank has polarity reversal
    for shank in channel_df.shank.unique():
        if channel_df[channel_df.shank == shank].channelClass.unique().shape[0] == 2:
            channel_df.loc[channel_df.shank == shank, "polarity_reversal"] = True
        else:
            channel_df.loc[channel_df.shank == shank, "polarity_reversal"] = False

    # add ripple and sharp wave features        
    labels = ["ripple_power", "ripple_amplitude", "SWR_diff", "SWR_amplitude"]
    for label in labels: 
        try:
            channel_df.loc[channel_sort_idx,label] = np.hstack(data[name][label][0][0][0])[0]
        except:
            x = np.arange(len(channel_sort_idx))*np.nan
            x[0:len(np.hstack(data[name][label][0][0][0])[0])] = np.hstack(data[name][label][0][0][0])[0]
            channel_df.loc[channel_sort_idx,label] = x

    # pull put avg ripple traces and ts
    ripple_time_axis = data[name]["ripple_time_axis"][0][0][0]
    ripple_average = np.ones([channel_df.shape[0],len(ripple_time_axis)])*np.nan

    rip_map = []
    for ch,values in zip(channels_,data[name]["ripple_average"][0][0][0]):
        if values.shape[1]>0:
            rip_map.append(values)
        else:
            rip_map.append(np.zeros([len(ripple_time_axis),len(ch)])*np.nan)
    
    ripple_average[channel_sort_idx] = np.hstack(rip_map).T

    brainRegions = load_brain_regions(basepath)
    for key, value in brainRegions.items():
        if 'ca1' in key.lower():
            for shank in value['electrodeGroups']:
                channel_df.loc[channel_df.shank == shank,"ca1_shank"] = True

    if (ripple_average.shape[0] != channel_df.shape[0]) & (~bypass_mismatch_exception):
        raise Exception('size mismatch '+
                        str(np.hstack(ripple_average).shape[1]) +
                        ' and ' +
                        str(channel_df.shape[0])
        )

    channel_df['basepath'] = basepath

    return channel_df, ripple_average, ripple_time_axis

def load_mua_events(basepath):
    """
    Loads the MUA data from the basepath.
    Meant to load .mat file created by find_HSE.m

    input:
        basepath: str
            The path to the folder containing the MUA data.
    output:
        mua_data: pandas.DataFrame
            The pandas.DataFrame containing the MUA data

    TODO: if none exist in basepath, create one
    """

    # locate .mat file
    try:
        filename = glob.glob(basepath+os.sep+'*mua_ca1_pyr.events.mat')[0]
    except:
        warnings.warn("file does not exist")
        return pd.DataFrame()

    # load matfile
    data = sio.loadmat(filename)

    # pull out and package data
    df = pd.DataFrame()
    df["start"] = data['HSE']["timestamps"][0][0][:, 0]
    df["stop"] = data['HSE']["timestamps"][0][0][:, 1]
    df["peaks"] = data['HSE']["peaks"][0][0]
    df["center"] = data['HSE']["center"][0][0]
    df["duration"] = data['HSE']["duration"][0][0]
    df["amplitude"] = data['HSE']["amplitudes"][0][0]
    df["amplitudeUnits"] = data['HSE']["amplitudeUnits"][0][0][0]
    df["detectorName"] = data['HSE']["detectorinfo"][0][0]["detectorname"][0][0][0]

    # get basename and animal
    normalized_path = os.path.normpath(filename)
    path_components = normalized_path.split(os.sep)
    df["basepath"] = basepath
    df["basename"] = path_components[-2]
    df["animal"] = path_components[-3]

    return df

def load_manipulation(
    basepath, struct_name=None, return_epoch_array=True, merge_gap=None
):
    """
    Loads the data from the basename.eventName.manipulations.mat file and returns a pandas dataframe

    file structure defined here:
        https://cellexplorer.org/datastructure/data-structure-and-format/#manipulations

    inputs:
        basepath: string, path to the basename.eventName.manipulations.mat file
        struct_name: string, name of the structure in the mat file to load. If None, loads all the manipulation files.
        return_epoch_array: bool, if True, returns only the epoch array
        merge_gap: int, if not None, merges the epochs that are separated by less than merge_gap (sec). return_epoch_array must be True.
    outputs:
        df: pandas dataframe, with the following columns:
            - start (float): start time of the manipulation in frames
            - stop (float): stop time of the manipulation in frames
            - peaks (float): list of the peak times of the manipulation in frames
            - center (float): center time of the manipulation in frames
            - duration (float): duration of the manipulation in frames
            - amplitude (float): amplitude of the manipulation
            - amplitudeUnit (string): unit of the amplitude
    Example:
        >> basepath = r"Z:\Data\Can\OML22\day8"
        >> df_manipulation = load_manipulation(basepath,struct_name="optoStim",return_epoch_array=False)
        >> df_manipulation.head(2)

                start	    stop	    peaks	    center	    duration amplitude amplitudeUnits
        0	8426.83650	8426.84845	8426.842475	8426.842475	0.01195	19651	pulse_respect_baseline
        1	8426.85245	8426.86745	8426.859950	8426.859950	0.01500	17516	pulse_respect_baseline

        >> basepath = r"Z:\Data\Can\OML22\day8"
        >> df_manipulation = load_manipulation(basepath,struct_name="optoStim",return_epoch_array=True)
        >> df_manipulation

        <EpochArray at 0x1faba577520: 5,774 epochs> of length 1:25:656 minutes
    """
    try:
        if struct_name is None:
            filename = glob.glob(basepath + os.sep + "*manipulation.mat")
            print(filename)
            if len(filename) > 1:
                raise ValueError(
                    "multi-file not implemented yet...than one manipulation file found"
                )
            filename = filename[0]
        else:
            filename = glob.glob(
                basepath + os.sep + "*" + struct_name + ".manipulation.mat"
            )[0]
    except:
        return None
    # load matfile
    data = sio.loadmat(filename)

    if struct_name is None:
        struct_name = list(data.keys())[-1]

    df = pd.DataFrame()
    df["start"] = data[struct_name]["timestamps"][0][0][:, 0]
    df["stop"] = data[struct_name]["timestamps"][0][0][:, 1]
    df["peaks"] = data[struct_name]["peaks"][0][0]
    df["center"] = data[struct_name]["center"][0][0]
    df["duration"] = data[struct_name]["duration"][0][0]
    df["amplitude"] = data[struct_name]["amplitude"][0][0]
    df["amplitudeUnits"] = data[struct_name]["amplitudeUnits"][0][0][0]

    if return_epoch_array:
        # get session epochs to add support for epochs
        epoch_df = load_epoch(basepath)
        # get session bounds to provide support
        session_bounds = nel.EpochArray(
            [epoch_df.startTime.iloc[0], epoch_df.stopTime.iloc[-1]]
        )
        manipulation_epoch = nel.EpochArray(
            np.array([df["start"], df["stop"]]).T, domain=session_bounds
        )
        if merge_gap is not None:
            manipulation_epoch = manipulation_epoch.merge(gap=merge_gap)

        return manipulation_epoch
    else:
        return df