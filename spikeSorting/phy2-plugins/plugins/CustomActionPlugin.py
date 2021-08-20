# import from plugins/action_status_bar.py
"""Show how to create new actions in the GUI.

The first action just displays a message in the status bar.

The second action selects the first N clusters, where N is a parameter that is entered by
the user in a prompt dialog.

"""

from phy import IPlugin, connect
import numpy as np
import logging

logger = logging.getLogger('phy')

class CustomActionPlugin(IPlugin):
    def attach_to_controller(self, controller):
        @connect
        def on_gui_ready(sender, gui):

            @controller.supervisor.actions.add(shortcut='ctrl+c')
            def select_first_unsorted():

                # All cluster view methods are called with a callback function because of the
                # asynchronous nature of Python-Javascript interactions in Qt5.
                @controller.supervisor.cluster_view.get_ids
                def find_unsorted(cluster_ids):
                    """This function is called when the ordered list of cluster ids is returned
                    by the Javascript view."""
                    groups = controller.supervisor.cluster_meta.get('group',list(range(max(cluster_ids))))
                    for ii in cluster_ids:
                        if groups[ii] == None or groups[ii] == 'unsorted':
                            s = controller.supervisor.clustering.spikes_in_clusters([ii])
                            if len(s)>0:
                                firstclu = ii
                                break
                    
                    if 'firstclu' in locals():
                        controller.supervisor.select(firstclu)

                    return


            @controller.supervisor.actions.add(shortcut='ctrl+v')
            def move_selected_to_end():

                logger.warn("Moving cluster to end")
                selected = controller.supervisor.selected
                s = controller.supervisor.clustering.spikes_in_clusters(selected)
                outliers2 = np.ones(len(s),dtype=int)
                controller.supervisor.actions.split(s,outliers2)
