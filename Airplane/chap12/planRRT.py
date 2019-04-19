import numpy as np
from message_types.msg_waypoints import msg_waypoints


class planRRT():
    def __init__(self):
        self.waypoints = msg_waypoints()
        self.segmentLength = 100 # standard length of path segments
        self.Distance = 5

    def planPath(self, wpp_start, wpp_end, map):

        # desired down position is down position of end node
        pd = wpp_end.item(2)

        # specify start and end nodes from wpp_start and wpp_end
        # format: N, E, D, cost, parentIndex, connectsToGoalFlag,
        start_node = np.array([[wpp_start.item(0), wpp_start.item(1), pd, 0, 0, 0]])
        end_node = np.array([[wpp_end.item(0), wpp_end.item(1), pd, 0, 0, 0]])

        # establish tree starting with the start node
        tree = start_node

        # check to see if start_node connects directly to end_node
        if ((np.linalg.norm(start_node[0:3] - end_node[0:3]) < self.segmentLength ) and not self.collision(start_node, end_node, map)):
            self.waypoints.ned = end_node[0:3]
        else:
            numPaths = 0
            while numPaths < 3:
                tree, flag = self.extendTree(tree, end_node, self.segmentLength, map, pd)
                numPaths = numPaths + flag


        # find path with minimum cost to end_node
        path = self.findMinimumPath(tree, end_node)
        # path = self.smoothPath(path,map)

        self.waypoints.type = 'fillet'
        self.waypoints.num_waypoints = path.shape[0]
        Va = 25
        self.waypoints.ned[:, 0:self.waypoints.num_waypoints] \
            = path[:,:3].T
        self.waypoints.airspeed[:, 0:self.waypoints.num_waypoints] \
            = Va*np.ones(self.waypoints.num_waypoints)

        return self.waypoints

    def generateRandomNode(self, map):
        n = np.random.uniform(low=0,high=map.city_width)
        e = np.random.uniform(low=0,high=map.city_width)
        return [n,e]

    def pointsAlongPath(self, start_node, end_node, Del):
        #Determine direction and distance between two points
        q = (end_node[0,:2] - start_node[0,:2]) / np.linalg.norm(end_node[0,:2] - start_node[0,:2])
        #Create array of points between the two points in that direction
        N = np.array([start_node.item(0)])
        E = np.array([start_node.item(1)])
        D = np.array([start_node.item(2)])
        for spot in range(1,int(np.round(np.linalg.norm(end_node[0,:2] - start_node[0,:2]),0)/self.Distance)):
            dis = q+self.Distance*spot
            N = np.append(N, start_node.item(0) + dis.item(0))
            E = np.append(E, start_node.item(1) + dis.item(1))
            D = np.append(D, start_node.item(2))
        N = np.append(N,end_node.item(0))
        E = np.append(E,end_node.item(1))
        D = np.append(D,start_node.item(2))
        #Return points
        return N, E, D

    def collision(self, start_node, end_node, map):
        radius = map.building_width/2.+5
        [N, E, D] = self.pointsAlongPath(start_node,end_node,self.Distance)
        for building_n in range(0, int(map.building_north.shape[0])):
            north_clearance = N-(map.building_north[building_n]-map.building_width/2.)
            north_bool = abs(north_clearance)<=radius
            H = -D[north_bool]
            if north_bool.any():
                # print("possible collision at:",building_n)
                for building_e in range(0, int(map.building_east.shape[0])):
                    east_clearance = E[north_bool] - (map.building_east[building_e] - map.building_width / 2.)
                    east_bool = abs(east_clearance) <= radius
                    if east_bool.any():
                        # print("possible collision at:", building_e)
                        H = H[east_bool]
                        H_clearance = H-map.building_height[building_n][building_e]
                        H_bool = H_clearance<=0
                        if H_bool.any():
                            # print("collision!!")
                            return True
        # print("No collisions")
        return False

        #
        # for spot in range(0,int(self.segmentLength/self.Distance)+1):
        #     # print("Checking for collision")
        #     #Check north point compared to north point pattern
        #
        #         if (abs(N[spot]-(map.building_north[building_n]-map.building_width/2.))<=radius):
        #             #Check east point compared to east point pattern
        #             for building_e in range(0, int(map.building_east.shape[0])):
        #                 if(abs(E[spot]-(map.building_east[building_e]-map.building_width/2.))<=radius):
        #                     # print("very possible collision")
        #                     if (-D[spot]<=map.building_height[building_n][building_e]):
        #                         # print("collision at spot:",spot,"!!")
        #                         return True
        #                 # print("east is safe at building:",building_e)
        #         # print("north is safe at building:",building_n)
        #     # print("segment spot:",spot," is safe.")
        # # print("No collisions on path")
        # return False

    def extendTree(self, tree, end_node, segmentLength, map, pd):
        # print("extendTree")
        NoCollision = False
        while NoCollision == False:
            #Generate random node
            n, e = self.generateRandomNode(map)
            # print("n:",n,"\ne:",e)
            #Find closest node already on the tree
            Dis = (n-tree[:,0])**2+(e-tree[:,1])**2
            new_node = np.array([[n,e]])
            spot = np.argmin(Dis)
            close_node = np.reshape(tree[spot,:],(1,6))
            #Create a node one segment length away from tree node
            q = (new_node[0,:2] - close_node[0,:2]) / np.linalg.norm(new_node[0,:2] - close_node[0,:2])
            new_node = close_node[0,:2]+self.segmentLength*q
            new_node = np.reshape(new_node,(1,2))
            Dis = np.linalg.norm(new_node[0,:2] - close_node[0,:2])
            # print("node:",new_node)
            if not self.collision(close_node,new_node,map):#If no collision
                NoCollision = True
                #Add node to tree with appropriate data
                cost = Dis+close_node[0,3]
                new_node = np.append(new_node,[end_node[0,2],cost,spot,0])
                new_node = new_node.reshape([1,6])
                tree = np.append(tree,new_node,axis=0)
                if (np.linalg.norm(tree[-1][0:3] - end_node[0,0:3]) < self.segmentLength):#if new node is one segment length away from final node
                    flag = 1
                    tree[-1][5]=1

                else:
                    flag = 0
        return [tree, flag]

    def findMinimumPath(self, tree, end_node):
        print("find Min Path")
        connection_spots = tree[:,-1]==1
        connections = tree[connection_spots]
        spot = np.argmin(connections[:,3])
        path = np.array([connections[spot]])
        spot = int(path[-1,4])
        prev_node = tree[spot,:]
        prev_node = np.reshape(prev_node,[1,6])
        path = np.append(path,prev_node,axis=0)
        while spot!=0:
            spot = int(prev_node[-1,4])
            prev_node = tree[spot,:]
            prev_node = np.reshape(prev_node,[1,6])
            path = np.append(path,prev_node,axis=0)
        path = np.flipud(path)
        return path

    def smoothPath(self, path, map):
        print(path.shape)
        Smoothpath = np.array([path[0,:]])
        current_node = np.reshape(Smoothpath[0, :],(1,6))
        prev_node = np.reshape(path[0, :],(1,6))
        for spot in range(1,np.shape(path)[0]):
            next_node = np.reshape(path[spot,:],(1,6))
            if self.collision(current_node,next_node,map):
                print("Updating smooth path at spot",spot)
                prev_node = np.reshape(prev_node,(1,6))
                Smoothpath = np.append(Smoothpath,prev_node,axis=0)
                current_node = prev_node
                prev_node = next_node
            else:
                prev_node = next_node
        print(Smoothpath)
        print("smoothPath")
        return path

