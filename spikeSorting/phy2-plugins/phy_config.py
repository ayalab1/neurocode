
# You can also put your plugins in ~/.phy/plugins/.

from phy import IPlugin

try:
    import phycontrib
except:
    pass

# Plugin example:
#
# class MyPlugin(IPlugin):
#     def attach_to_cli(self, cli):
#         # you can create phy subcommands here with click
#         pass

c = get_config()
c.Plugins.dirs = [r'~/.phy/plugins/']
c.TemplateGUI.plugins = ['Recluster','SplitShortISI'] 
# ['AmplitudeHistogram','ControllerSettings','Recluster','Export2neurosuite'] # c.KwikGUI.plugins = ['ControllerSettingsKwik']