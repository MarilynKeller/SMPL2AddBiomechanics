from aitviewer.utils.colors import vertex_colors_from_weights
import numpy as np

# osim_markers_color = ()
osim_mesh_color = (160 / 255, 160 / 255, 160 / 255, 1.0)

# osso_markers_color = ()
osso_mesh_color = (150/255, 150/255, 50/255, 1)

mocap_markers_color = ()

skin_mesh_color = (179/255, 120/255, 179/255, 0.5)
skin_mesh_oppaque_color = (149/255, 85/255, 149/255, 0.9)

skel_color = (60/255,80/255,200/255,1)


def color_per_index(mesh):
    """ Color a trimesh mesh according to its vertex indices."""
    mesh.visual.vertex_colors[:] = vertex_colors_from_weights(np.arange(len(mesh.vertices)), shuffle=False, alpha=1)*255
    return mesh 