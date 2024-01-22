# Copyright (C) 2024  Max Planck Institute for Intelligent Systems Tuebingen, Marilyn Keller 

import numpy as np
import trimesh

def compute_mapping_barycentric(verts, mesh):
    """ 
    For each vertex, compute the closest point on the mesh and return:
    - the closest point on the mesh
    - the closest triangle on the mesh
    - the coordinates of the vertex in this triangle frame

    vertices: numpy array of shape (n_vertices, 3)
    mesh: trimesh.base.Trimesh object
    output: mesh_verts_idx, mesh_tri_index, verts_coords
    """

    #Find closest triangle on osso
    proximityQuery = trimesh.proximity.ProximityQuery(mesh)
    _, _, mesh_tri_index  = proximityQuery.on_surface(verts) 

    # import ipdb; ipdb.set_trace()
    triangles = mesh.triangles[mesh_tri_index]
    verts_coords = trimesh.triangles.points_to_barycentric(triangles, verts, method='cramer') 

    return mesh_tri_index, verts_coords

def apply_mapping_barycentric(baricentric_verts_coords, mesh_tri_index, mesh):
    """
    Apply the mapping computed by compute_mapping to the mesh
    """
    import ipdb; ipdb.set_trace()
    triangles = mesh.triangles[mesh_tri_index]
    normals, _ = trimesh.triangles.normals(triangles)
    verts_coords =  trimesh.triangles.barycentric_to_points(triangles, baricentric_verts_coords[:,:]) # + baricentric_verts_coords[:,-1, np.newaxis] * normals
    return verts_coords
 


def compute_mapping(verts, mesh):
    """ 
    For each vertex, compute the closest point on the mesh and return:
    - the closest point on the mesh
    - the closest triangle on the mesh
    - the coordinates of the vertex in this triangle frame

    vertices: numpy array of shape (n_vertices, 3)
    mesh: trimesh.base.Trimesh object
    output: mesh_verts_idx, mesh_tri_index, verts_coords
    """

    #Find closest triangle on osso
    proximityQuery = trimesh.proximity.ProximityQuery(mesh)
    _, _, mesh_tri_index  = proximityQuery.on_surface(verts) 


    # We want to represent the point m_vect_abs[i] (x,y,z) in the base attached to the closest triangle (x',y',z,)
    # Note that this base is not orthogonal
    # v0 = self.osso_p.vertices[self.markers_osso_idx]

    # Find triangle points 

    # First we build the base (u,v,w) from the triangle
    triangles = mesh.triangles[mesh_tri_index]
    v0 = triangles[:,0]
    v1 = triangles[:,1]
    v2 = triangles[:,2]
    
    m_vect_abs = verts - v0 #skeleton to marker vectors in posed space

    u = (v1-v0)/np.linalg.norm(v1-v0, axis=1)[:,np.newaxis]
    v = (v2-v0)/np.linalg.norm(v2-v0, axis=1)[:,np.newaxis]
    w = np.cross(u, v) / np.linalg.norm(np.cross(u, v), axis=1)[:,np.newaxis]

    # Then we compute the passage matrix such that (x,y,z) = A.T * (x',y',z')
    # (https://public.iutenligne.net/mathematiques/algebre/arrou-vignod/changement_de_base/exprimer_un_vecteur_v_dans_la_nouvelle_base.html)

    A = np.dstack((u, v, w)) # Nx3x3

    # I did not manage to do the matrix operation per i in (0:N) for I do it in a loop
    m_vect_rel = np.zeros_like(m_vect_abs) # skeleton to marker vectors in triangle space
    for i in range(m_vect_abs.shape[0]):
        m_vect_rel[i] = np.matmul(A[i].T, m_vect_abs[i])

    # TEST
    # Test reconstruction, the result should be the same as m_vect_abs
    # The equation used is inv(A.T) * (x,y,z) = (x',y',z')
    m_rec_test = np.zeros_like(m_vect_abs)          
    for i in range(m_vect_abs.shape[0]):
        m_rec_test[i] = np.matmul(np.linalg.inv(A[i].T), m_vect_rel[i])

    m_rec_test - m_vect_abs
    assert(np.linalg.norm(m_rec_test - m_vect_abs)< 1e-10), f'Error in mapping: {np.linalg.norm(m_rec_test - m_vect_abs)}'

    return mesh_tri_index, m_vect_rel



def apply_mapping(mesh_tri_index, m_vect_rel, mesh):
    """
    Apply the mapping to recover markers from the mesh
    """
    # Recover triangles
    triangles = mesh.triangles[mesh_tri_index]
    v0 = triangles[:,0]
    v1 = triangles[:,1]
    v2 = triangles[:,2]

    # Build the per triangle frame of reference
    u = (v1-v0)/np.linalg.norm(v1-v0, axis=1)[:,np.newaxis]
    v = (v2-v0)/np.linalg.norm(v2-v0, axis=1)[:,np.newaxis]
    w = np.cross(u, v) / np.linalg.norm(np.cross(u, v), axis=1)[:,np.newaxis]

    # Then we compute the passage matrix such that (x,y,z) = A.T * (x',y',z')
    # (https://public.iutenligne.net/mathematiques/algebre/arrou-vignod/changement_de_base/exprimer_un_vecteur_v_dans_la_nouvelle_base.html)

    A = np.dstack((u, v, w)) # Nx3x3

    # transform the mesh to marker vector from triangle frame to world frame
    m_vect_abs = np.zeros_like(m_vect_rel)          
    for i in range(m_vect_rel.shape[0]):
       m_vect_abs[i] = np.matmul(np.linalg.inv(A[i].T), m_vect_rel[i])

    m_reconstructed = v0  + m_vect_abs

    return m_reconstructed


def get_submesh(verts, faces, verts_retained=None, faces_retained=None, min_vert_in_face=2):
    '''
        Given a mesh, create a (smaller) submesh
        indicate faces or verts to retain as indices or boolean

        @return new_verts: the new array of 3D vertices
                new_faces: the new array of faces
                bool_faces: the faces indices wrt the input mesh
                vetex_ids: the vertex_ids wrt the input mesh
        '''

    # import ipdb; ipdb.set_trace()
    if verts_retained is not None:
        # Transform indices into bool array
        if verts_retained.dtype != 'bool':
            vert_mask = np.zeros(len(verts), dtype=bool)
            vert_mask[verts_retained] = True
        else:
            vert_mask = verts_retained

        # Faces with at least min_vert_in_face vertices
        bool_faces = np.sum(vert_mask[faces.ravel()].reshape(-1, 3), axis=1) > min_vert_in_face

    elif faces_retained is not None:
        # Transform indices into bool array
        if faces_retained.dtype != 'bool':
            bool_faces = np.zeros(len(faces_retained), dtype=bool)
        else:
            bool_faces = faces_retained

    new_faces = faces[bool_faces]
    # just in case additional vertices are added
    vertex_ids = np.where(vert_mask)[0]

    oldtonew = -1 * np.ones([len(verts)])
    oldtonew[vertex_ids] = range(0, len(vertex_ids))

    new_verts = verts[vertex_ids]
    new_faces = oldtonew[new_faces].astype('int32')

    # if new_verts.shape[0] != len(verts_retained):
    #     import ipdb; ipdb.set_trace()

    return (new_verts, new_faces, bool_faces, vertex_ids)

