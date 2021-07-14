#read data from txt, and output the file
class Read_txt:
    def __init__(self,based_profit):
        self.based_profit=based_profit
        self.file="1.txt"
        self.physical_depot_id=0
        self.g_number_of_vehicles=5
        self.vehicle_arrival_time_ending=300
        print("Reading data...")

    def read_nodes(self):
        self.node_list=[]
        self.g_number_of_nodes=0
        self.g_number_of_commodities=0
        self.compartment_capacity=0
        with open(self.file,"r") as fl:
            lines=fl.readlines()
            self.g_number_of_commodities=int(lines[1])
            self.g_number_of_nodes = int((int(lines[0])-1)/self.g_number_of_commodities+1)
            list = lines[2].strip().split(" ")
            self.compartment_capacity=int(list[0])
            index = 1
            for line in lines[4:4+int(lines[0])]:
                list=line.strip().split(" ")
                if int(list[0])<self.g_number_of_nodes:
                    node=Node()
                    node.node_id=int(list[0])
                    node.x=float(list[1])
                    node.y=float(list[2])
                    node.demand=[int(list[3])]
                    node.base_profit_for_lagrangian=[self.based_profit]*self.g_number_of_commodities
                    node.base_profit_for_searching=[self.based_profit]*self.g_number_of_commodities
                    self.node_list.append(node)
                else:
                    # node_id=int(list[0])%(self.g_number_of_nodes-1)#node_id
                    node=self.node_list[index]
                    node.demand.append(int(list[3]))
                    index+=1
                    if index==self.g_number_of_nodes:
                        index=1
            # print()

    def generate_links(self):
        self.read_nodes()
        self.link_list=[]
        self.g_number_of_links=0
        for from_node in range(0,self.g_number_of_nodes):
            for to_node in range(0, self.g_number_of_nodes):
                if from_node == to_node:
                    continue
                else:
                    x_from_node = self.node_list[from_node].x
                    y_from_node = self.node_list[from_node].y
                    x_to_node = self.node_list[to_node].x
                    y_to_node = self.node_list[to_node].y
                    distance=((x_from_node-x_to_node)**2+(y_from_node-y_to_node)**2)**0.5

                    link=Link()
                    link.link_id=self.g_number_of_links
                    link.from_node_id=from_node
                    link.to_node_id=to_node

                    link.distance=distance
                    self.link_list.append(link)

                    self.node_list[from_node].outbound_nodes_list.append(to_node)
                    self.node_list[from_node].outbound_links_list.append(link)
                    self.node_list[from_node].outbound_size=len(self.node_list[from_node].outbound_links_list)

                    self.g_number_of_links+=1

        node = Node()
        node.node_id =self.g_number_of_nodes
        self.virtural_depot =self.g_number_of_nodes # 增加虚拟车场
        node.x = 0
        node.y = 0
        node.demand = 0
        node.base_profit_for_lagrangian =None
        node.base_profit_for_searching = None
        self.node_list.append(node)
        self.g_number_of_nodes+=1

        link = Link() # 增加路段 （虚拟-实际车场）
        link.link_id = self.g_number_of_links

        link.from_node_id = self.virtural_depot
        link.to_node_id = self.physical_depot_id

        link.distance =1
        self.link_list.append(link)

        self.node_list[link.from_node_id].outbound_nodes_list.append(link.to_node_id)
        self.node_list[link.from_node_id].outbound_links_list.append(link)
        self.node_list[link.from_node_id].outbound_size = len(self.node_list[link.from_node_id].outbound_links_list)

        self.g_number_of_links += 1

        link = Link() # 增加路段 （实际-虚拟车场）
        link.link_id = self.g_number_of_links

        link.from_node_id = self.physical_depot_id
        link.to_node_id = self.virtural_depot

        link.distance = 1
        self.link_list.append(link)

        self.node_list[link.from_node_id].outbound_nodes_list.append(link.to_node_id)
        self.node_list[link.from_node_id].outbound_links_list.append(link)
        self.node_list[link.from_node_id].outbound_size = len(self.node_list[link.from_node_id].outbound_links_list)

        self.g_number_of_links += 1
        # self.generate_GAMS_file()
        return self.node_list,self.link_list,self.g_number_of_nodes,self.g_number_of_commodities,self.compartment_capacity\
            ,self.g_number_of_vehicles,self.vehicle_arrival_time_ending,self.virtural_depot,self.physical_depot_id

    def generate_GAMS_file(self):
        with open("GAMS_input.txt", "w") as fl:
            # index
            fl.write("set i node id /0*{}/;\n".format(self.g_number_of_nodes))
            fl.write("set v vehicle id /1*{}/;\n".format(self.g_number_of_vehicles))
            fl.write("set t time id /0*{}/;\n".format(self.vehicle_arrival_time_ending))
            fl.write("set f fraction id /0*{}/;\n".format(self.g_number_of_commodities))
            fl.write("alias(i,j);\n")
            fl.write("alias(t,s);\n")

            # parameters-arc cost
            fl.write("parameter arcs(i,j,t,s) arc cost ;\n")
            for link in self.link_list:
                arc_cost = link.distance
                from_node = "'" + str(link.from_node_id) + "',"
                to_node = "'" + str(link.to_node_id) + "',"
                fl.write("arcs(" + from_node + to_node + "t,t+{})={};\n".format(arc_cost, arc_cost))

            node = "'" + str(self.physical_depot_id) + "',"
            fl.write("arcs(" + node + node + "t,t+{})={};\n".format(1, 0.0001))

            #demand
            fl.write("Parameter demand(i,f) demand ;\n")
            for node_id in range(1,self.g_number_of_nodes):
                node=self.node_list[node_id]
                demand_list=node.demand
                for c in range(self.g_number_of_commodities):
                    demand=demand_list[c]
                    node_str="'"+str(node_id)+"',"
                    f_str="'"+str(c)+"'"
                    fl.write("demand("+node_str+f_str+")={};\n".format(demand))

            # parameters-capacity
            fl.write("parameter capacity(v,f);\n")
            fl.write("capacity(v,f)={};\n".format(self.compartment_capacity))

            # flow balance

            fl.write("parameter origin_node(v,i,t)  origin nodes and departure time; \n")
            fl.write("origin_node(v,'{}','{}') = 1;\n".format(self.physical_depot_id, 0))

            fl.write("parameter destination_node(v,i,t);\n")
            fl.write("destination_node(v,'{}','{}') = 1;\n".format(self.physical_depot_id, self.vehicle_arrival_time_ending))

            fl.write("parameter intermediate_node(v,i,t);\n")
            fl.write("intermediate_node(v,i,t) = (1- origin_node(v,i,t))*(1- destination_node(v,i,t));\n")
            print("OK!")

    def output_file(self):
        with open("Node.csv","w") as fl:
            fl.write("Node_id,x_coordinate,y_coordinate,demand\n")
            for node in self.node_list:
                node_id=node.node_id
                x_coordinate=node.x
                y_coordinate=node.y
                demand=node.demand
                fl.write(str(node_id)+","+str(x_coordinate)+","+str(y_coordinate)+","+str(demand)+"\n")

        with open("link.csv","w") as fl:
            fl.write("link_id,from_node,to_node,distance\n")
            for link in self.link_list:
                link_id=link.link_id
                from_node=link.from_node_id
                to_node=link.to_node_id
                distance=link.distance
                fl.write(str(link_id)+","+str(from_node)+","+str(to_node)+","+str(distance)+"\n")

class Node:
    def __init__(self):
        self.node_id=None
        self.node_type=None
        self.x=None
        self.y=None
        self.demand=None
        self.outbound_nodes_list=[]
        self.outbound_size=0
        self.outbound_links_list=[]
        self.base_profit_for_searching = None  # 三项之力 rpo*η-rpo/2+Π（pai）
        self.base_profit_for_lagrangian = None  # 乘子 Π（pai）

class Link:
    def __init__(self):
        self.link_id=None
        # self.link_type=None     #物理路段-1、虚拟路段-0
        self.from_node_id=None
        self.to_node_id=None
        self.distance=None

# mod=Read_txt(0)
# mod.generate_links()
# mod.output_file()