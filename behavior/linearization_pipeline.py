import glob
import os
import pickle
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from track_linearization import make_track_graph
from track_linearization import get_linearized_position
from scipy.io import savemat, loadmat

plt.ion()
plt.style.use("dark_background")

"""
TODO: 
-fix bugs in removing nodes and edges
-automatically determine get_linearized_position params "sensor_std_dev"
    for instances where data is not in cm
-add option to input your own nodes and edges
"""


class NodePicker:
    """Interactive creation of track graph by looking at video frames."""

    def __init__(
        self, ax=None, basepath=None, node_color="#177ee6", node_size=100, epoch=None
    ):
        if ax is None:
            ax = plt.gca()
        self.ax = ax
        self.canvas = ax.get_figure().canvas
        self.cid = None
        self._nodes = []
        self.node_color = node_color
        self._nodes_plot = ax.scatter([], [], zorder=5, s=node_size, color=node_color)
        self.edges = [[]]
        self.basepath = basepath
        self.epoch = epoch

        if self.epoch is not None:
            self.epoch = int(self.epoch)

        ax.set_title(
            "Left click to place node.\nRight click to remove node."
            "\nShift+Left click to clear nodes.\nCntrl+Left click two nodes to place an edge"
            "\nEnter to save and exit.",
            fontsize=8,
        )

        self.canvas.draw()

        self.connect()

    @property
    def node_positions(self):
        return np.asarray(self._nodes)

    def connect(self):
        if self.cid is None:
            self.cid = self.canvas.mpl_connect("button_press_event", self.click_event)
            self.canvas.mpl_connect("key_press_event", self.process_key)

    def disconnect(self):
        if self.cid is not None:
            self.canvas.mpl_disconnect(self.cid)
            self.cid = None

    def process_key(self, event):
        if event.key == "enter":
            self.format_and_save()

    def click_event(self, event):
        if not event.inaxes:
            return
        if (event.key not in ["control", "shift"]) & (event.button == 1):  # left click
            self._nodes.append((event.xdata, event.ydata))
        if (event.key not in ["control", "shift"]) & (event.button == 3):  # right click
            self.remove_point((event.xdata, event.ydata))
        if (event.key == "shift") & (event.button == 1):
            self.clear()
        if (event.key == "control") & (event.button == 1):
            point = (event.xdata, event.ydata)
            distance_to_nodes = np.linalg.norm(self.node_positions - point, axis=1)
            closest_node_ind = np.argmin(distance_to_nodes)
            if len(self.edges[-1]) < 2:
                self.edges[-1].append(closest_node_ind)
            else:
                self.edges.append([closest_node_ind])
        if event.key == "enter":
            self.format_and_save()

        self.redraw()

    def redraw(self):
        # Draw Node Circles
        if len(self.node_positions) > 0:
            self._nodes_plot.set_offsets(self.node_positions)
        else:
            self._nodes_plot.set_offsets([])

        # Draw Node Numbers
        for ind, (x, y) in enumerate(self.node_positions):
            self.ax.text(
                x,
                y,
                ind,
                zorder=6,
                fontsize=10,
                horizontalalignment="center",
                verticalalignment="center",
                clip_on=True,
                bbox=None,
                transform=self.ax.transData,
            )
        # Draw Edges
        for edge in self.edges:
            if len(edge) > 1:
                x1, y1 = self.node_positions[edge[0]]
                x2, y2 = self.node_positions[edge[1]]
                self.ax.plot(
                    [x1, x2], [y1, y2], color="#1f8e4f", linewidth=3, zorder=1000
                )
        self.canvas.draw()

    def remove_point(self, point):
        if len(self._nodes) > 0:
            distance_to_nodes = np.linalg.norm(self.node_positions - point, axis=1)
            closest_node_ind = np.argmin(distance_to_nodes)
            self._nodes.pop(closest_node_ind)

    def clear(self):
        self._nodes = []
        self.edges = [[]]
        self.redraw()

    def format_and_save(self):

        behave_df = load_animal_behavior(self.basepath)

        if self.epoch is not None:
            epochs = load_epoch(self.basepath)
            # na_idx = np.isnan(behave_df.x) | (
            #     (behave_df.time < epochs.iloc[self.epoch].startTime)
            #     & (behave_df.time > epochs.iloc[self.epoch].stopTime)
            # )
            cur_epoch = (
                ~np.isnan(behave_df.x) &
                (behave_df.time >= epochs.iloc[self.epoch].startTime) &
                (behave_df.time <= epochs.iloc[self.epoch].stopTime)
            )
        else:
            # na_idx = np.isnan(behave_df.x)
            cur_epoch = ~np.isnan(behave_df.x)

        print("running hmm...")
        track_graph = make_track_graph(self.node_positions, self.edges)

        position = np.vstack(
            [behave_df[cur_epoch].x.values, behave_df[cur_epoch].y.values]
        ).T

        position_df = get_linearized_position(
            position=position,
            track_graph=track_graph,
            edge_order=self.edges,
            use_HMM=True,
        )

        print("saving to disk...")
        behave_df.loc[cur_epoch, "linearized"] = position_df.linear_position.values
        behave_df.loc[cur_epoch, "states"] = position_df.track_segment_id.values
        behave_df.loc[cur_epoch, "projected_x_position"] = position_df.projected_x_position.values
        behave_df.loc[cur_epoch, "projected_y_position"] = position_df.projected_y_position.values

        filename = glob.glob(os.path.join(self.basepath, "*.animal.behavior.mat"))[0]
        data = loadmat(filename, simplify_cells=True)

        data["behavior"]["position"]["linearized"] = behave_df.linearized.values
        data["behavior"]["states"] = behave_df.states.values
        data["behavior"]["position"]["projected_x"] = behave_df.projected_x_position.values
        data["behavior"]["position"]["projected_y"] = behave_df.projected_y_position.values

        # store nodes and edges within behavior file
        data = self.save_nodes_edges_to_behavior(data, behave_df)

        savemat(filename, data, long_field_names=True)

        self.save_nodes_edges()
        self.disconnect()
        plt.close()

    def save_nodes_edges(self):
        results = {"node_positions": self.node_positions, "edges": self.edges}
        save_file = os.path.join(self.basepath, "linearization_nodes_edges.pkl")
        with open(save_file, "wb") as f:
            pickle.dump(results, f)

    def save_nodes_edges_to_behavior(self, data, behave_df):
        """
        Store nodes and edges into behavior file
        Searches to find epochs with valid linearized coords
        Nodes and edges are stored within behavior.epochs{n}.{node_positions and edges}
        """ 
        if self.epoch is None:

            # load epochs
            epochs = load_epoch(self.basepath)
            # iter over each epoch
            for epoch_i, ep in enumerate(epochs.itertuples()):
                # locate index for given epoch
                idx = behave_df.time.between(ep.startTime, ep.stopTime)
                # if linearized is not all nan, add nodes and edges
                if not all(np.isnan(behave_df[idx].linearized)) & (behave_df[idx].shape[0] != 0):
                    # adding nodes and edges
                    data["behavior"]["epochs"][epoch_i]["node_positions"] = self.node_positions
                    data["behavior"]["epochs"][epoch_i]["edges"] = self.edges
        else:
            # if epoch was used, add nodes and edges just that that epoch
            data["behavior"]["epochs"][self.epoch]["node_positions"] = self.node_positions
            data["behavior"]["epochs"][self.epoch]["edges"] = self.edges

        return data
        
def load_animal_behavior(basepath):
    filename = glob.glob(os.path.join(basepath, "*.animal.behavior.mat"))[0]
    data = loadmat(filename, simplify_cells=True)
    df = pd.DataFrame()
    df["time"] = data["behavior"]["timestamps"]
    try:
        df["states"] = data["behavior"]["states"]
    except:
        pass
    for key in data["behavior"]["position"].keys():
        try:
            df[key] = data["behavior"]["position"][key]
        except:
            pass
    return df


def load_epoch(basepath):
    """
    Loads epoch info from cell explorer basename.session and stores in df
    """
    filename = glob.glob(os.path.join(basepath, "*.session.mat"))[0]

    data = loadmat(filename, simplify_cells=True)
    try:
        return pd.DataFrame(data["session"]["epochs"])
    except:
        return pd.DataFrame([data["session"]["epochs"]])


def run(basepath, epoch=None):
    print("here is the file,", basepath)
    fig, ax = plt.subplots(figsize=(5, 5))

    behave_df = load_animal_behavior(basepath)

    if epoch is not None:
        epochs = load_epoch(basepath)

        behave_df = behave_df[
            behave_df["time"].between(
                epochs.iloc[epoch].startTime, epochs.iloc[epoch].stopTime
            )
        ]

    ax.scatter(behave_df.x, behave_df.y, color="white", s=0.5, alpha=0.5)
    ax.axis("equal")
    ax.set_axisbelow(True)
    ax.yaxis.grid(color="gray", linestyle="dashed")
    ax.xaxis.grid(color="gray", linestyle="dashed")
    ax.set_ylabel("y (cm)")
    ax.set_xlabel("x (cm)")

    picker = NodePicker(ax=ax, basepath=basepath, epoch=epoch)

    plt.show(block=True)


if __name__ == "__main__":
    print(len(sys.argv))
    if len(sys.argv) == 2:
        run(sys.argv[1])
    elif len(sys.argv) == 3:
        run(sys.argv[1], epoch=int(sys.argv[2]))