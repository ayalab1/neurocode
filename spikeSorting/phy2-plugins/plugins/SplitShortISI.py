"""Show how to write a custom split action."""

from phy import IPlugin, connect
import numpy as np
import logging

logger = logging.getLogger('phy')



class SplitShortISI(IPlugin):
    def attach_to_controller(self, controller):
        @connect
        def on_gui_ready(sender, gui):
            #@gui.edit_actions.add(shortcut='alt+i')
            @controller.supervisor.actions.add(shortcut='alt+i')
            def VisualizeShortISI():
                """Split all spikes with an interspike interval of less than 1.5 ms into a separate
                cluster. THIS IS FOR VISUALIZATION ONLY, it will show you where potential noise
                spikes may be located. Re-merge the clusters again afterwards and cut the cluster with
                another method!"""
                
                logger.info('Detecting spikes with ISI less than 1.5 ms')

                # Selected clusters across the cluster view and similarity view.
                cluster_ids = controller.supervisor.selected

                # Get the amplitudes, using the same controller method as what the amplitude view
                # is using.
                # Note that we need load_all=True to load all spikes from the selected clusters,
                # instead of just the selection of them chosen for display.
                bunchs = controller._amplitude_getter(cluster_ids, name='template', load_all=True)

                # We get the spike ids and the corresponding spike template amplitudes.
                # NOTE: in this example, we only consider the first selected cluster.
                spike_ids = bunchs[0].spike_ids
                spike_times = controller.model.spike_times[spike_ids]
                dspike_times = np.diff(spike_times)

                labels = np.ones(len(dspike_times),'int64')
                labels[dspike_times<.0015]=2
                labels = np.append(labels,1) #include last spike to match with len spike_ids

                # We perform the clustering algorithm, which returns an integer for each
                # subcluster.
                #labels = k_means(y.reshape((-1, 1)))
                assert spike_ids.shape == labels.shape

                # We split according to the labels.
                controller.supervisor.actions.split(spike_ids, labels)
                logger.info('Splitted short ISI spikes from main cluster')
