"""
Function for reading .2dm file and creating "points, vertices, boundary" output required to create an Anuga domain/mesh
-- arg1: mesh.2dm file path/name
-- arg2: boundaries in order of assignment in the mesh.2dm file, e.g.: 'rivright', 'ds', 'rivleft', 'us'

From ANUGA Manual (page 6), the items needed to build mesh from 'first principles are:'
-- a list "points" giving the coordinates of each mesh point,
-- a list "vertices" specifying the three vertices of each triangle, and
-- a dictionary "boundary" that stores the edges on the boundary and associates with each a symbolic tag. The edges are
   represented as pairs (i, j) where i refers to the triangle id and j to the edge id of that triangle. Edge ids are
   enumerated from 0 to 2 based on the id of the vertex opposite.

"""


def Convert_mesh(smsFile, bndys):

    with open(smsFile) as mesh2dm:
        #mesh2dm.readline()
        #mesh2dm.readline()
        lines = mesh2dm.readlines()

    points, vertices, boundary = [], [], {}  # these are outputs required by Anuga domain method
    boundaryNodes = {}
    for bndy in bndys:
        boundaryNodes[bndy] = []  # create dictionary of empty lists, each bndy tag will become a key
    tag = len(bndys)

    # iterate through all lines in the files (minus first two text lines)
    for line in lines:
        s = line.strip().split()

        # the .2dm file starts with vertex IDs defining each triangle
        # this reads those IDs into a list of tuples, where the "minus 1" is necessary to convert to python zero-indexing
        if s[0] == 'E3T':
            vertices.append((int(s[2])-1, int(s[3])-1, int(s[4])-1))

        # add vertex coordinate tuples (x_coord, y_coord) to list
        if s[0] == 'ND':
            points.append((float(s[2]), float(s[3])))

        # add boundary nodes to dictionary of lists where bndy tags are keys
        if s[0] == 'NS':
            bndy = bndys[len(bndys) - tag]
            for node in s[1:]:
                ni = int(node)
                if ni < 0:
                    ni = ni*(-1)  # last node in nodestring is indicated by negative ID
                    boundaryNodes[bndy].append(ni-1)
                    tag -= 1
                    # number following last node ID is the nodestring ID, this is not a node ID so we break to avoid it
                    break
                boundaryNodes[bndy].append(ni-1)  # subtract 1 to adjust to zero-indexing

    # this block creates a list of indices of elements that have two vertices on the boundary...
    # which means the element itself is on the boundary
    # we need to identify edges that are on the boundary (see definition above)
    index = 0  # index variables keeps track of element ID (b/c "vertices" list is ordered)
    boundaryNodesJoined = sum(boundaryNodes.values(), [])  # create list of just values (boundary vertex IDs), w/o tags
    for i, j, k in vertices:
        elemOnBndy = False

        # start with a check for corner elements, at the junction of two nodestrings, where all three verts are on boundary
        # if true, our element is at a corner, and has two edges on boundary
        # this is unusual; more of a just in case step
        if all(v in boundaryNodesJoined for v in [i, j, k]):  # need to look through ALL bndy nodes for this
            elemOnBndy, edgeID, tags = True, [], []
            for tag in bndys:  # need to find which TWO boundaries
                if all(v in boundaryNodes[tag] for v in [i, j]):
                    # apparently edge ID is enumerated from 0 to 2 based on the ID of the opposite vertex,
                    # i.e., the one not on the boundary
                    # so, if vert IDs i and j (0 and 1) are on bndy, vert k (2) is not, and thus the bndy edge ID is 2
                    edgeID.append(2)
                    tags.append(tag)
                if all(v in boundaryNodes[tag] for v in [j, k]):
                    edgeID.append(0)  # means node 0/i is opposite this boundary edge
                    tags.append(tag)
                if all(v in boundaryNodes[tag] for v in [k, i]):
                    edgeID.append(1)
                    tags.append(tag)

        # now we check for elements with one edge on the boundary (vast majority of cases)
        # if 2 of 3 nodes are on the boundary (as defined by nodestrings), the element is on the boundary
        else:
            for tag in bndys:  # loop through the boundary tags
                if all(v in boundaryNodes[tag] for v in [i, j]):
                    elemOnBndy, edgeID = True, 2
                    # stop looking once we find the tag we want, tag variable is applied to "boundary" dict at end
                    break
                elif all(v in boundaryNodes[tag] for v in [j, k]):  
                    elemOnBndy, edgeID = True, 0
                    break        
                elif all(v in boundaryNodes[tag] for v in [k, i]): 
                    elemOnBndy, edgeID = True, 1
                    break

        # add boundary dict entry, using edgeID + tag from above, and current element ID (index)
        if elemOnBndy:
            # add boundary tags passed from arg2 to each boundary element
            if type(edgeID) is list:
                boundary[(index, edgeID[0])] = tags[0]
                boundary[(index, edgeID[1])] = tags[1]
            else:
                boundary[(index, edgeID)] = tag
                
        index += 1

    return points, vertices, boundary
