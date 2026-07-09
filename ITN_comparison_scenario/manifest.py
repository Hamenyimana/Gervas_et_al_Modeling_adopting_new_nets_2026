import os
dirname = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sif_dir = os.path.dirname(__file__)

schema_file      = os.path.join(dirname, "download/schema.json")
eradication_path = os.path.join(dirname, "download/Eradication")
assets_input_dir = "Assets"
plugins_folder   = "download/reporter_plugins"
sif_path         = os.path.join(sif_dir, 'dtk_sif.id')
my_ep4_assets    = None
