from Read_data import Read_txt
import copy
import itertools
import math
class LR:
    def __init__(self):
        base_profit=-20
        data=Read_txt(base_profit)
        self.node_list,self.link_list,self.g_number_of_nodes,self.g_number_of_commodities,self.compartment_capacity\
            ,self.g_number_of_vehicles,self.vehicle_arrival_time_ending,self.virtural_depot,self.physical_depot_id=\
            data.generate_links()

        self.vehicle_departure_time_beginning = 0
        self.big_M=1000
        self.iteration_times=100
        self.rpo=1
        self.generate_node_serving_state_set()
        #beam-search DP
        self.g_ending_state_vector = [None]*(self.g_number_of_vehicles+1)#多加一辆车给上界用
        self.g_time_dependent_state_vector = []


        for v in range((self.g_number_of_vehicles+1)):
            self.g_time_dependent_state_vector.append([])
            intervals = self.vehicle_arrival_time_ending-self.vehicle_departure_time_beginning
            for t in range(intervals+1):
                # self.g_time_dependent_state_vector[v].append([None]*self.g_number_of_nodes)
                self.g_time_dependent_state_vector[v].append([None])

        # SST path 要输出的路径*4
        self.path_node_seq = []
        self.path_time_seq = []
        self.path_commodity_state_seq = []
        self.path_vehicle_state_seq = []

        self.record_profit = [] # 要输出的乘子
        # gap
        self.served_time = []
        self.repeat_served = []
        self.un_served=[]
        self.max_label_cost = float("inf")

        #Bounds
        self.ADMM_local_LB = [0] * self.iteration_times
        self.ADMM_local_UB = [0] * self.iteration_times
        self.ADMM_global_LB = [-self.max_label_cost] * self.iteration_times # 要输出的上下界+gap
        self.ADMM_global_UB = [self.max_label_cost] * self.iteration_times


    def g_solve_CSMC_by_LR(self):
        for i in range(self.iteration_times): #从（0）开始到（迭代次数-1）
            self.path_node_seq.append([])
            self.path_time_seq.append([])
            self.path_commodity_state_seq.append([])
            self.path_vehicle_state_seq.append([])
            self.record_profit.append([])

            serving_times=[]
            for node_id in range(self.g_number_of_nodes):
                serving_times.append([0]*self.g_number_of_commodities)

            self.served_time.append(serving_times)

            serving_times_1 = []
            for node_id in range(self.g_number_of_nodes):
                serving_times_1.append([0] * self.g_number_of_commodities)
            self.repeat_served.append(serving_times_1)

            serving_times_2 = []
            for node_id in range(self.g_number_of_nodes):
                serving_times_2.append([0] * self.g_number_of_commodities)
            self.un_served.append(serving_times_2)

            # self.record_profit.append(serving_times)
            print("iteration: {}".format(i+1))# 索引加1，才是从第1次开始迭代，不然输出的是0-39（迭代次数是40的话）
#step1:UB generation
            if i != 0:
                self.served_time[i] = copy.copy(self.served_time[i-1]) #复制上次迭代的服务次数信息
            for v in range(self.g_number_of_vehicles):
                if self.g_ending_state_vector[v]!= None:
                    path_commodity_state_seq=self.g_ending_state_vector[v].VSStateVector[0].path_commodity_state_seq
                    for node_id in range(1,self.g_number_of_nodes-1):#只对1-10这10个客户点，self.g_number_of_nodes-1=11，但右边是开区间取不到11
                        for c in range(self.g_number_of_commodities):
                            self.served_time[i][node_id][c]-=path_commodity_state_seq[node_id][c] # 要拿掉当前车的服务次数才是第i次迭代的η

                for node_id in range(1,self.g_number_of_nodes-1):
                    for c in range(self.g_number_of_commodities):
                        Node=self.node_list[node_id]
                        Node.base_profit_for_searching[c]=Node.base_profit_for_lagrangian[c]+self.rpo*self.served_time[i][node_id][c]-self.rpo/2 #对每辆车更新每个commodity的高级乘子

                self.g_optimal_time_dependent_dynamic_programming(v, 1) # 求解每辆车的上界解
                self.ADMM_local_UB[i]+=self.g_ending_state_vector[v].VSStateVector[0].Primal_Label_cost # 局部上界 目标函数 += cx

                self.path_node_seq[i].append(self.g_ending_state_vector[v].VSStateVector[0].path_node_seq)
                self.path_time_seq[i].append(self.g_ending_state_vector[v].VSStateVector[0].path_time_seq)
                self.path_commodity_state_seq[i].append(self.g_ending_state_vector[v].VSStateVector[0].path_commodity_state_seq)
                self.path_vehicle_state_seq[i].append(self.g_ending_state_vector[v].VSStateVector[0].path_vehicle_state_seq)

                path_commodity_state_seq = self.g_ending_state_vector[v].VSStateVector[0].path_commodity_state_seq
                for node_id in range(1,self.g_number_of_nodes-1):
                    for c in range(self.g_number_of_commodities):
                        self.served_time[i][node_id][c] += path_commodity_state_seq[node_id][c] # 之前拿掉了是为了计算η，在这里要再加回来，给下一次迭代用

            #Check
            Flag_check=0
            for node_id in range(1,self.g_number_of_nodes-1):
                for c in range(self.g_number_of_commodities):
                    if self.served_time[i][node_id][c]>1:  #重复访问
                        Flag_check=1                       #记录到repeat_served列表
                        self.repeat_served[i][node_id][c]=1
                    if self.served_time[i][node_id][c]==0: #没有访问
                        self.ADMM_local_UB[i]+=self.big_M  #惩罚
                        self.un_served[i][node_id][c] = 1  #记录到un_served列表
            # if Flag_check==0:
            #     print(self.ADMM_local_UB[i])
            print("serving times:")
            print(self.served_time[i])


# step2: LB generation
            self.g_optimal_time_dependent_dynamic_programming(self.g_number_of_vehicles, 2) #求所有车的下界解，每一辆车都一样
            if len(self.g_ending_state_vector[self.g_number_of_vehicles].VSStateVector)<self.g_number_of_vehicles:
                self.ADMM_local_LB[i]=self.g_number_of_vehicles*self.g_ending_state_vector[self.g_number_of_vehicles].VSStateVector[0].Label_cost_for_lagrangian #计算下界的第一项，子问题的目标函数
            else: # 最终集合里有多于车辆数条路径的话，按成本大小逐条分配，这样下界会高一点
                for v in range(self.g_number_of_vehicles):
                    self.ADMM_local_LB[i] += self.g_ending_state_vector[self.g_number_of_vehicles].VSStateVector[v].Label_cost_for_lagrangian
            #add multipliers
            for node_id in range(1,self.g_number_of_nodes-1):
                for c in range(self.g_number_of_commodities):
                    multiplier=self.node_list[node_id].base_profit_for_lagrangian[c] #拿出每个commodity的普通乘子
                    self.ADMM_local_LB[i]-=multiplier #加上下界的第二项，-sum(i)sum(f)pai

            # print(self.ADMM_local_LB[i])

# step3: multiplier updating
            if i==0:
                self.ADMM_global_LB[i]=self.ADMM_local_LB[i]
            else:
                self.ADMM_global_LB[i] = max(self.ADMM_local_LB[i],self.ADMM_global_LB[i-1])

            if i==0:
                if Flag_check==0: # 全cover1次或含未服务的
                    self.ADMM_global_UB[i] = self.ADMM_local_UB[i]
                else:
                    self.ADMM_global_UB[i]=10000
            else: #第二次迭代开始
                if Flag_check == 0:# 全cover1次或含未服务的
                    self.ADMM_global_UB[i] = min(self.ADMM_local_UB[i],self.ADMM_global_UB[i-1])
                else:# Flag_check=1，即有重复服务的
                    self.ADMM_global_UB[i] =self.ADMM_global_UB[i-1]

            for node_id in range(1,self.g_number_of_nodes-1):
                multiplier = []
                for c in range(self.g_number_of_commodities):
                    Node = self.node_list[node_id]
                    multiplier.append(Node.base_profit_for_lagrangian[c])
                    Node.base_profit_for_lagrangian[c]+=self.rpo*(self.served_time[i][node_id][c]-1) #
                self.record_profit[i].append(multiplier)# 记录乘子

            # 第i次迭代
            unserved_num_i=0
            repeated_num_i=0
            for node_id in range(1,self.g_number_of_nodes-1):
                for c in range(self.g_number_of_commodities):
                    if self.un_served[i][node_id][c]==1:
                        unserved_num_i+=1
                    if self.repeat_served[i][node_id][c]==1:
                        repeated_num_i+=1
            # 第i-1次迭代
            unserved_num_i_1 = 0
            repeated_num_i_1= 0
            for node_id in range(1, self.g_number_of_nodes - 1):
                for c in range(self.g_number_of_commodities):
                    if self.un_served[i-1][node_id][c] == 1:
                        unserved_num_i_1 += 1
                    if self.repeat_served[i-1][node_id][c] == 1:
                        repeated_num_i_1 += 1
            # 第i次迭代，惩罚参数的更新方法
            if i >= 10:
                if (unserved_num_i+repeated_num_i)**2>0.25 * (unserved_num_i_1+repeated_num_i_1) ** 2:#第i次的解不好
                    self.rpo += 1 #调大乘子
                if (unserved_num_i+repeated_num_i)**2 == 0:#第i次的解还行，没有重复访问，也没有未访问的
                    self.rpo = 0.5


    def g_optimal_time_dependent_dynamic_programming(self,vehicle_id,Flag):

        self.origin_node=self.virtural_depot
        self.destination_node=self.virtural_depot
        for t in range(self.vehicle_departure_time_beginning, self.vehicle_arrival_time_ending+1):
            # for n in range(self.g_number_of_nodes):
            self.g_time_dependent_state_vector[vehicle_id][t] = C_time_indexed_state_vector()
            self.g_time_dependent_state_vector[vehicle_id][t].Reset()
            self.g_time_dependent_state_vector[vehicle_id][t].current_time=t
            # self.g_time_dependent_state_vector[vehicle_id][t].current_node=n

        self.g_ending_state_vector[vehicle_id]=C_time_indexed_state_vector()


        element = CVSState(self.g_number_of_nodes,self.g_number_of_commodities)
        element.current_node_id=self.origin_node
        element.current_time=self.vehicle_departure_time_beginning
        element.current_commodity_serving_state=[0]*self.g_number_of_commodities
        element.current_vehicle_state=[self.compartment_capacity]*self.g_number_of_commodities

        element.path_node_seq.append(self.origin_node)
        element.path_time_seq.append(self.vehicle_departure_time_beginning)
        element.path_vehicle_state_seq.append([self.compartment_capacity]*self.g_number_of_commodities)
        # element.path_commodity_state_seq.append([0]*self.g_number_of_commodities)
        self.g_time_dependent_state_vector[vehicle_id][self.vehicle_departure_time_beginning].update_state(element,Flag)
        # DP
        for t in range(self.vehicle_departure_time_beginning, self.vehicle_arrival_time_ending):
            self.g_time_dependent_state_vector[vehicle_id][t].Sort(Flag)
            for index in range(min(len(self.g_time_dependent_state_vector[vehicle_id][t].VSStateVector),100)):
                pElement=self.g_time_dependent_state_vector[vehicle_id][t].VSStateVector[index] #pElement是CVSState的实例化，C_time_indexed_state_vector？
                from_node_id=pElement.current_node_id
                from_node=self.node_list[from_node_id]
                for i in range(from_node.outbound_size):
                    to_node_id=from_node.outbound_nodes_list[i]
                    to_node=self.node_list[to_node_id]
                    link_to=from_node.outbound_links_list[i]
                    next_time=t+math.ceil(link_to.distance)
                    if next_time>self.vehicle_arrival_time_ending:
                        continue
                    # if pElement.path_node_serving_state[to_node_id]==1:
                    #     continue

                    # Case 1: destination
                    if to_node_id==self.destination_node:
                        new_element=CVSState(self.g_number_of_nodes,self.g_number_of_commodities)
                        new_element.my_copy(pElement)

                        new_element.current_node_id=self.destination_node
                        new_element.current_time=next_time
                        new_element.current_commodity_serving_state=[0]*self.g_number_of_commodities
                        # new_element.current_vehicle_state=[self.compartment_capacity]*self.g_number_of_commodities

                        new_element.path_node_seq.append(self.destination_node)
                        new_element.path_time_seq.append(next_time)
                        new_element.path_vehicle_state_seq.append(new_element.current_vehicle_state)
                        # new_element.path_commodity_state_seq.append([0]*self.g_number_of_commodities)

                        state=[0]*self.g_number_of_commodities
                        # new_element.Calculate_Label_Cost(to_node, link_to, state)

                        # using the update can reduce the scale of VSState
                        self.g_ending_state_vector[vehicle_id].VSStateVector.append(new_element)
                        continue

                    # Case 2:physical depot
                    if to_node_id == self.physical_depot_id:
                        new_element = CVSState(self.g_number_of_nodes, self.g_number_of_commodities)
                        new_element.my_copy(pElement)

                        new_element.current_node_id = self.physical_depot_id
                        new_element.current_time = next_time

                        new_element.current_commodity_serving_state = [0] * self.g_number_of_commodities
                        # new_element.current_vehicle_state=[self.compartment_capacity]*self.g_number_of_commodities

                        new_element.path_node_seq.append(self.physical_depot_id)
                        new_element.path_time_seq.append(next_time)
                        new_element.path_vehicle_state_seq.append(new_element.current_vehicle_state)
                        # new_element.path_commodity_state_seq.append([0]*self.g_number_of_commodities)

                        state = [0] * self.g_number_of_commodities
                        if from_node_id != self.origin_node:
                            new_element.Calculate_Label_Cost(to_node, link_to, state)

                        self.g_time_dependent_state_vector[vehicle_id][next_time].update_state(new_element, Flag)
                        continue


                    # Case 3: customer node
                    if pElement.path_node_serving_state[to_node_id]==1:
                        continue
                    for state in self.state_list:
                        #check capacity
                        Flag_=0
                        residual_ca_list=[]
                        for c in range(self.g_number_of_commodities):
                            residual_ca=pElement.current_vehicle_state[c]-to_node.demand[c]*state[c]
                            residual_ca_list.append(residual_ca)
                            if state[c]==1 and residual_ca<0:
                                Flag_=1
                        if Flag_==1:
                            continue

                        #transition is ok
                        new_element=CVSState(self.g_number_of_nodes,self.g_number_of_commodities)
                        new_element.my_copy(pElement)

                        new_element.current_node_id=to_node_id
                        new_element.current_time=next_time
                        new_element.current_commodity_serving_state=state

                        new_element.current_vehicle_state=residual_ca_list

                        new_element.path_node_seq.append(to_node_id)
                        new_element.path_time_seq.append(next_time)
                        new_element.path_vehicle_state_seq.append(new_element.current_vehicle_state)

                        # new_element.path_commodity_state_seq.append(state)
                        new_element.path_commodity_state_seq[to_node_id]=state

                        new_element.path_node_serving_state[to_node_id]=1

                        new_element.Calculate_Label_Cost(to_node,link_to,state)
                        self.g_time_dependent_state_vector[vehicle_id][next_time].update_state(new_element,Flag)
                        continue

        self.g_ending_state_vector[vehicle_id].Sort(Flag)

        print(self.g_ending_state_vector[vehicle_id].VSStateVector[0].path_node_seq)
        print(self.g_ending_state_vector[vehicle_id].VSStateVector[0].path_commodity_state_seq)


    def generate_node_serving_state_set(self):
        self.state_list= [list(i) for i in itertools.product([0, 1], repeat=self.g_number_of_commodities)]
        # print(self.state_list)
        self.state_list.pop(0)
        # print(self.state_list)

    def output_to_file(self, spend_time):
        # GAP
        with open("Gap.csv", "w") as f:
            f.write("iteration,LB,UB,gap,served_time\n")
            for i in range(self.iteration_times):
                f.write(str(i + 1) + ",")
                f.write(str(self.ADMM_global_LB[i]) + ",")
                f.write(str(self.ADMM_global_UB[i]) + ",")
                gap=(self.ADMM_global_UB[i]-self.ADMM_global_LB[i])/self.ADMM_global_UB[i]
                f.write(str(round(gap,3)) + ",")

                served_time=self.served_time[i]
                str1=""
                for node_id in range(self.g_number_of_nodes):
                    str2=""
                    for serving_time in served_time[node_id]:
                        str2+=str(serving_time)+"+"
                    str1+=str2[:-1]+"_"
                f.write(str1[:-1])
                f.write("\n")
            f.write("time: {} seconds".format(round(spend_time,2)))

        # 乘子
        with open("Multipliers.csv", "w") as f:
            f.write("Iteration" + ",")
            for node_id in range(1,self.g_number_of_nodes-1):
                for c in range(self.g_number_of_commodities):
                    f.write(str(node_id)+"_"+str(c+1) + ",")
            f.write("\n")

            for i in range(self.iteration_times):
                f.write(str(i+1) + ",")
                for node_id in range(0, self.g_number_of_nodes - 2):
                    for c in range(self.g_number_of_commodities):
                        multiplier=self.record_profit[i][node_id][c]
                        f.write(str(multiplier)+",")
                f.write("\n")


        # SST路径
        with open("Path.csv", "w") as f:
            f.write("iteration,vehicle_id,path_node_seq,path_time_seq,path_vehicle_state_seq,path_commodity_state_seq\n")
            for i in range(self.iteration_times):
                for v in range(self.g_number_of_vehicles):
                    path_node_seq=self.path_node_seq[i][v]
                    path_time_seq=self.path_time_seq[i][v]
                    path_commodity_state_seq=self.path_commodity_state_seq[i][v]
                    path_vehicle_state_seq=self.path_vehicle_state_seq[i][v]

                    commodity_state=[]
                    for node_id in path_node_seq:
                        commodity_state.append(path_commodity_state_seq[node_id])

                    f.write(str(i + 1) + ",")
                    f.write(str(v + 1) + ",")
                    str1 = ""
                    str2 = ""
                    str3 = ""
                    str4 = ""

                    for s in range(len(path_node_seq)):
                        state_str=""
                        for state in commodity_state[s]:
                            state_str+=str(state)+"+" #commodity的state

                        state_str_1 = ""
                        for state in path_vehicle_state_seq[s]:
                            state_str_1 += str(state) + "+" #compartment的state

                        str1 += str(path_node_seq[s]) + "_"
                        str2 += str(path_time_seq[s]) + "_"
                        str3 += state_str_1[:-1] + "_" # vehicle_state
                        str4+=state_str[:-1]+"_" # commodity_state
                    f.write(str1[:-1] + "," + str2[:-1] + "," + str3[:-1] +","+str4[:-1]+"\n")


class C_time_indexed_state_vector:
    def __init__(self):
        self.current_time=0
        # self.current_node=0
        self.VSStateVector=[]
        self.state_map=[] #node+[node serving state]

    def Reset(self):
        self.current_time =0
        # self.current_node =0
        self.VSStateVector=[]
        self.state_map=[]

    def m_find_state_index(self, string_key):  #tring_key为水量
        if string_key in self.state_map:
            return self.state_map.index(string_key)
        else:
            return -1

    def update_state(self, element, Flag):
        string_key = element.generate_string_key()
        state_index = self.m_find_state_index(string_key)
        if state_index==-1:#True:
            self.VSStateVector.append(element)
            self.state_map.append(string_key)
        else:
            # Flag=1,ADMM成本；Flag=2,LR成本
            if Flag==2:
                if element.Label_cost_for_lagrangian<self.VSStateVector[state_index].Label_cost_for_lagrangian:
                    self.VSStateVector[state_index]=element

            if Flag==1:
                if element.Label_cost_for_searching<self.VSStateVector[state_index].Label_cost_for_searching:
                    self.VSStateVector[state_index]=element


    def Sort(self, Flag):
        if Flag == 1:
            self.VSStateVector = sorted(self.VSStateVector, key=lambda x: x.Label_cost_for_searching)
        if Flag == 2:
            self.VSStateVector = sorted(self.VSStateVector, key=lambda x: x.Label_cost_for_lagrangian)


class CVSState:
    def __init__(self, g_number_of_nodes,g_number_of_commodity):
        #node state
        self.current_node_id=0
        self.current_time=0

        self.current_commodity_serving_state=[]
        self.current_vehicle_state=[]

        #path state
        self.path_node_seq=[]
        self.path_time_seq=[]


        self.path_commodity_state_seq = []
        for node_id in range(g_number_of_nodes):
            self.path_commodity_state_seq.append([0]*g_number_of_commodity)

        self.path_vehicle_state_seq=[]

        self.path_node_serving_state=[0]*g_number_of_nodes

        self.Primal_Label_cost=0
        self.Label_cost_for_lagrangian=0
        self.Label_cost_for_searching=0

    def generate_string_key(self):#"1_[1,0]"
        str_=""
        str_+=str(self.current_node_id)+"-"
        str_+=str(self.current_commodity_serving_state)
        return str_

    def my_copy(self, pElement):
        #node state
        self.current_node_id=copy.copy(pElement.current_node_id)
        self.current_time=copy.copy(pElement.current_time)

        self.current_commodity_serving_state=[]
        self.current_commodity_serving_state=copy.copy(pElement.current_commodity_serving_state)

        self.current_vehicle_state=[]
        self.current_vehicle_state=copy.copy(pElement.current_vehicle_state)

        #path state
        self.path_node_seq=[]
        self.path_node_seq=copy.copy(pElement.path_node_seq)

        self.path_time_seq=[]
        self.path_time_seq=copy.copy(pElement.path_time_seq)

        self.path_commodity_state_seq = []
        self.path_commodity_state_seq=copy.copy(pElement.path_commodity_state_seq)

        self.path_vehicle_state_seq=[]
        self.path_vehicle_state_seq=copy.copy(pElement.path_vehicle_state_seq)

        self.path_node_serving_state=[]
        self.path_node_serving_state=copy.copy(pElement.path_node_serving_state)

        self.Primal_Label_cost=copy.copy(pElement.Primal_Label_cost)
        self.Label_cost_for_lagrangian=copy.copy(pElement.Label_cost_for_lagrangian)
        self.Label_cost_for_searching=copy.copy(pElement.Label_cost_for_searching)

    def Calculate_Label_Cost(self,node,link,state):

        self.Primal_Label_cost+=link.distance
        self.Label_cost_for_lagrangian += link.distance
        self.Label_cost_for_searching += link.distance

        for c in range(len(state)):
            if state[c]==1:
                self.Label_cost_for_lagrangian+=node.base_profit_for_lagrangian[c]
                self.Label_cost_for_searching+=node.base_profit_for_searching[c]